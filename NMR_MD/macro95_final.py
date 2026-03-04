# -*- coding: utf-8 -*-
"""
CCPNmr AnalysisAssign v3
NMR–MD GUI — (NMRbox)
- Input tab:
- Preprocessing tab:
- Minimisation/Relaxation tab adds:
  * Option to include NMR distance restraints (DISANG) if an RST file is present
  * CPU count for pmemd.MPI stages (heating/relaxation)
- Production tab:
- Analysis tab

"""

import datetime
import json
import os
import pathlib
import re
import shlex
import shutil
import stat
import subprocess

from ccpn.ui.gui.popups.Dialog import CcpnDialog
from ccpn.ui.gui.widgets.Button import Button
from ccpn.ui.gui.widgets.CheckBox import CheckBox
from ccpn.ui.gui.widgets.Entry import Entry, FloatEntry
from ccpn.ui.gui.widgets.Frame import Frame
from ccpn.ui.gui.widgets.HLine import HLine
from ccpn.ui.gui.widgets.Label import Label
from ccpn.ui.gui.widgets.MessageDialog import showWarning
from ccpn.ui.gui.widgets.PulldownList import PulldownList
from ccpn.ui.gui.widgets.RadioButtons import RadioButtons
from PyQt5.QtCore import QProcess
from PyQt5.QtWidgets import QFileDialog, QPlainTextEdit, QTabWidget

SHELL = "/bin/bash"
DEFAULT_VENV = "/home/nmrbox/0063/hkolli/nmrmd_venv"
AMBER_SH = "/usr/software/amber/amber.sh"

PREP_PDB_SH = r"""#!/usr/bin/env bash
set -euo pipefail

# Default placeholders (avoid unbound under set -u)
NP="${NP:-16}"
NPROCS="${NPROCS:-16}"

# Self-contained AMBER MD pre-processing pipeline (GUI-driven)
# Implements the steps from Hima's prep_pdb.sh, but accepts GUI-style flags.
#
# Flags:
#   --workdir <dir>                 working directory (expects source/input.pdb, optional source/input.nef)
#   --ph <float>                    pH for pdb2pqr/propka
#   --salt <float>                  salt molarity (M)
#   --protein_ff <ff>               e.g. ff19SB, ff14SB (leaprc.protein.<ff>)
#   --water <water>                 e.g. tip3p, spce, opc (leaprc.water.<water>)
#   --box <float>                   buffer distance in Å
#   --shape <oct|cube>              oct or cube
#
# Notes:
# - NEF conversion is optional and runs only if source/input.nef exists.
# - Neutralisation is derived from net charge computed from the protonated PDB.
# - Salt ion pairs computed from tleap-reported volume (Å^3).

PH="7.0"
SALT_M="0.150"
PROT_FF="ff19SB"
WATER="tip3p"
CUSHION_A="12.0"
SHAPE="oct"
WORKDIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --workdir) WORKDIR="$2"; shift 2;;
    --ph) PH="$2"; shift 2;;
    --salt) SALT_M="$2"; shift 2;;
    --protein_ff) PROT_FF="$2"; shift 2;;
    --water) WATER="$2"; shift 2;;
    --box) CUSHION_A="$2"; shift 2;;
    --shape) SHAPE="$2"; shift 2;;
    *) echo "WARNING: ignoring unknown arg '$1'"; shift;;
  esac
done


[[ -z "$WORKDIR" ]] && { echo "ERROR: --workdir required"; exit 2; }

SRC_DIR="${WORKDIR}/source"
INPDB="${SRC_DIR}/input.pdb"
NEF="${SRC_DIR}/input.nef"

[[ -f "$INPDB" ]] || { echo "ERROR: missing ${INPDB}. Save input first."; exit 2; }

OUT="${WORKDIR}/01_preprocess"
mkdir -p "$OUT"
# Choose preprocess staging dir without creating nested 01_preprocess/01_preprocess
if [ "$OUT" = "01_preprocess" ]; then
  PREP_DIR="$OUT"
else
  PREP_DIR="$OUT/01_preprocess"
fi
mkdir -p "$PREP_DIR"
# Stage inputs into OUT so relative paths in minrelax.sh work when running inside OUT
for f in input_salt.prmtop input_salt.inpcrd input.RST; do
  if [ -f "./01_preprocess/$f" ]; then
    cp -f "./01_preprocess/$f" "$PREP_DIR/$f"
  elif [ -f "./$f" ]; then
    cp -f "./$f" "$PREP_DIR/$f"
  fi
done
# Copy helper script into OUT
if [ -f "./minrelax.sh" ]; then
  cp -f "./minrelax.sh" "$OUT/minrelax.sh"
  chmod +x "$OUT/minrelax.sh"
fi

# Ensure NPROCS comes from Condor allocation
# NPROCS: prefer 1st arg from HTCondor (we pass $(request_cpus)), else env, else nproc.
if [ -n "${1:-}" ] && [[ "${1}" =~ ^[0-9]+$ ]]; then
  NPROCS="${1}"
elif [ -n "${NPROCS:-}" ]; then
  NPROCS="${NPROCS}"
elif [ -n "${_CONDOR_SLOT_CPUS:-}" ]; then
  NPROCS="${_CONDOR_SLOT_CPUS}"
else
  NPROCS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || nproc 2>/dev/null || echo 1)"
fi
echo "NPROCS=$NPROCS" >> "$OUT/_job_meta.txt"
cd "$OUT"
ROOT="$(basename "$INPDB" .pdb)"   # "input" by default

UNIT="SYS"
ION_POS="Na+"
ION_NEG="Cl-"

FF_PROT="leaprc.protein.${PROT_FF}"
WATER_LEAP="leaprc.water.${WATER}"

# water box name mapping for tleap solvent boxes
case "${WATER}" in
  tip3p) WATBOX="TIP3PBOX" ;;
  spce|spc) WATBOX="SPCBOX" ;;
  opc) WATBOX="OPCBOX" ;;
  *) WATBOX="${WATER^^}BOX" ;;
esac

run_cmd(){ echo -e "\n>>> $*"; "$@"; }

# --- Check dependencies ---
for cmd in pdb4amber tleap cpptraj python3; do
  command -v "$cmd" >/dev/null || { echo "ERROR: $cmd not found"; exit 1; }
done
if command -v pdb2pqr30 >/dev/null; then P2PQR=pdb2pqr30
elif command -v pdb2pqr >/dev/null; then P2PQR=pdb2pqr
else echo "ERROR: pdb2pqr (or pdb2pqr30) not found"; exit 1; fi

# --- Embedded Python functions ---
net_charge_py(){ python3 - "$1" <<'PY'
import sys
from collections import defaultdict

AA=set("ALA ARG ASN ASP CYS GLN GLU GLY HIS HID HIE HIP ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
RC={"ASP":-1,"GLU":-1,"LYS":1,"ARG":1,"HID":0,"HIE":0,"HIP":1,"HIS":0}
CAP={"ACE","NME"}; WAT={"HOH","WAT","SOL"}
ION={"NA":1,"K":1,"LI":1,"CL":-1,"BR":-1,"I":-1,"MG":2,"CA":2,"ZN":2,"MN":2,"CO":2,"NI":2,"CU":2,"CD":2,"SR":2,"BA":2,"FE":2,"FE2":2,"FE3":3,"AL":3,"YB":3,"HG":2}

def is_w(r):
  r=r.upper()
  return r in WAT or r.startswith("TIP")

def parse(p):
  R,Ch={},defaultdict(list)
  with open(p) as fh:
    for L in fh:
      if not (L.startswith("ATOM") or L.startswith("HETATM")):
        continue
      a=L[12:16].strip()
      r=L[17:20].strip().upper()
      c=(L[21].strip() or "_")
      i=L[22:26].strip()
      ic=L[26].strip()
      k=(c,i,ic)
      if k not in R:
        R[k]={"res":r,"atoms":set(),"is_pro":r in AA,"is_w":is_w(r),"is_cap":r in CAP,"is_ion":r in ION}
        Ch[c].append(k)
      R[k]["atoms"].add(a)
  for c in Ch:
    Ch[c]=sorted(Ch[c], key=lambda k:(c, int(k[1]) if k[1].isdigit() else 9999, k[2]))
  return R,Ch

def his_q(atoms):
  s={x.upper() for x in atoms}
  return 1 if ("HD1" in s and "HE2" in s) else 0

def charge(R,Ch):
  t=0
  for k,v in R.items():
    rn=v["res"]
    if v["is_w"]:
      continue
    if v["is_ion"]:
      t += ION.get(rn,0)
      continue
    if v["is_cap"]:
      continue
    if v["is_pro"]:
      if rn in RC:
        t += RC[rn]
      elif rn=="HIS":
        t += his_q(v["atoms"])
  # termini (very simple heuristic)
  for c,ks in Ch.items():
    n=[R[k]["res"] for k in ks]
    p=[i for i,k in enumerate(ks) if R[k]["is_pro"]]
    if not p:
      continue
    f,l=p[0],p[-1]
    if not (f-1>=0 and n[f-1]=="ACE"):
      t += 1
    if not (l+1<len(ks) and n[l+1]=="NME"):
      t -= 1
  return t

R,Ch=parse(sys.argv[1])
print(int(charge(R,Ch)))
PY
}

ion_pairs_py(){ python3 - "$1" "$2" <<'PY'
import sys
AVOG=6.02214076e23
A3L=1e-27
V=float(sys.argv[1])
C=float(sys.argv[2])
print(int(round(C*V*A3L*AVOG)))
PY
}

extract_volume(){
  awk '/[Vv]olume[[:space:]]*[:=]/ {
         for(i=1;i<=NF;i++) if($i ~ /^[0-9.]+$/) v=$i
       } END{if(v!="") print v}' "$1"
}

# 0) Optional NEF -> RST/DIP (runs only if NEF exists)
if [[ -f "$NEF" ]]; then
  if command -v nef_to_RST >/dev/null 2>&1; then
    echo "==> Step 0: Convert NEF → RST/DIP (will run after pdb4amber cleanup)"
  else
    echo "WARNING: nef_to_RST not found; skipping NEF conversion."
  fi
else
  echo "==> No NEF file provided — skipping NEF to RST conversion"
fi

# 1) pdb4amber cleanup
echo "==> Step 1: pdb4amber cleanup"
run_cmd pdb4amber -i "$INPDB" -o "${ROOT}_new.pdb" > pdb4amber_input.log 2>&1
echo "    Wrote ${ROOT}_new.pdb"

# 2) remove CRYST1 (robust)
echo "==> Step 2: Removing CRYST1 line(s)"
awk 'sub(/\r$/,""){}; toupper(substr($0,1,6))!="CRYST1"' "${ROOT}_new.pdb" > "${ROOT}_nocrys.pdb"
if grep -n '^[[:space:]]*CRYST1' "${ROOT}_nocrys.pdb" >/dev/null 2>&1; then
  echo "ERROR: CRYST1 still present after filtering. Aborting." >&2
  exit 1
fi

# 0) NEF conversion now that we have a cleaned pdb
# Detect pdb4amber renumber-map file (name depends on pdb4amber version)
P4A_MAP=""
if ls -1 *_renum.txt >/dev/null 2>&1; then
  P4A_MAP="$(ls -1 *_renum.txt | head -n1)"
elif ls -1 *renum*.txt >/dev/null 2>&1; then
  P4A_MAP="$(ls -1 *renum*.txt | head -n1)"
fi

if [[ -f "$NEF" ]] && command -v nef_to_RST >/dev/null 2>&1; then
  echo "==> Step 0 (continued): NEF → RST/DIP"
  # Create an AMBER-normalised PDB for nef_to_RST (more robust than using cleaned PDB directly)
  cat > leap_nef.in <<EOF
set default PBRadii mbondi3
source ${FF_PROT}
x = loadpdb ${ROOT}_new.pdb
saveAmberParm x ${ROOT}_nef.parm7 ${ROOT}_nef.rst7
quit
EOF
  tleap -f leap_nef.in > tleap_nef.out 2>&1
  ambpdb -p ${ROOT}_nef.parm7 -c ${ROOT}_nef.rst7 > ${ROOT}.ambpdb.pdb

  set +e
  if [[ -n "$P4A_MAP" && -f "$P4A_MAP" ]]; then
    nef_to_RST -nef "$NEF" -pdb "${ROOT}.ambpdb.pdb" -p4a "$P4A_MAP" -rst "${ROOT}.RST" -rdc "${ROOT}.DIP" > nef_to_RST.log 2>&1
  else
    nef_to_RST -nef "$NEF" -pdb "${ROOT}.ambpdb.pdb" -rst "${ROOT}.RST" -rdc "${ROOT}.DIP" > nef_to_RST.log 2>&1
  fi
  rc=$?
  set -e

  if [[ $rc -ne 0 ]]; then
    echo "WARNING: nef_to_RST exited with code $rc. See nef_to_RST.log"
  else
    if [[ ! -s "${ROOT}.RST" ]]; then
      echo "WARNING: ${ROOT}.RST is empty; restraints will be disabled. Check nef_to_RST.log"
    else
      echo "    Wrote ${ROOT}.RST and ${ROOT}.DIP"
    fi
  fi
fi

# 3) protonate using pdb2pqr + propka
echo "==> Step 3: PDB2PQR protonation (pH=${PH})"
run_cmd "$P2PQR" "${ROOT}_nocrys.pdb" "${ROOT}.pqr" \
  --pdb-output "${ROOT}_pqr.pdb" \
  --ffout AMBER --with-ph="${PH}" --titration-state-method=propka \
  > pdb2pqr.log 2>&1
echo "    Wrote ${ROOT}_pqr.pdb"

# 4) compute net charge (from protonated PDB)
echo "==> Step 4: Compute net charge"
NET_Q="$(net_charge_py "${ROOT}_pqr.pdb")"; [[ -z "$NET_Q" ]] && NET_Q=0
ABS_Q=$(python3 -c "print(abs(int('$NET_Q')))")
if (( NET_Q < 0 )); then NEUT_ION="$ION_POS"; NEUT_COUNT="$ABS_Q"
elif (( NET_Q > 0 )); then NEUT_ION="$ION_NEG"; NEUT_COUNT="$ABS_Q"
else NEUT_ION="$ION_POS"; NEUT_COUNT=0; fi
printf "    Net charge %+d → Neutralize with %d %s\n" "$NET_Q" "$NEUT_COUNT" "$NEUT_ION"

# 5) solvate + neutralize (tleap) and parse volume
echo "==> Step 5: LEaP solvation + neutralisation"
SOLVATE_CMD="solvatebox"
if [[ "$SHAPE" == "oct" ]]; then
  SOLVATE_CMD="solvateOct"
fi

cat > leap_solvate.in <<EOF
set default PBRadii mbondi3
source ${FF_PROT}
source ${WATER_LEAP}
${UNIT}=loadpdb ${ROOT}_pqr.pdb
addIons ${UNIT} ${NEUT_ION} ${NEUT_COUNT}
${SOLVATE_CMD} ${UNIT} ${WATBOX} ${CUSHION_A}
saveAmberParm ${UNIT} ${ROOT}_neutral.prmtop ${ROOT}_neutral.inpcrd
quit
EOF

run_cmd tleap -f leap_solvate.in | tee leap_solvate.log

V_A3="$(extract_volume leap_solvate.log)"
if [[ -z "$V_A3" ]]; then
  sides=$(awk '/Total[[:space:]]+bounding[[:space:]]+box[[:space:]]+for[[:space:]]+atom[[:space:]]+centers:/{
    for(i=1;i<=NF;i++) if($i~/^[0-9.]+$/){printf("%s ",$i)}}' leap_solvate.log | tail -n1)
  if [[ -n "$sides" ]]; then
    V_A3=$(python3 - "$sides" <<'PY'
import sys
a,b,c=map(float,sys.argv[1].split())
print(f"{a*b*c:.6f}")
PY
)
  fi
fi

[[ -z "$V_A3" ]] && { echo "ERROR: Could not parse volume from tleap log"; exit 1; }
echo "    Parsed tleap volume: ${V_A3} Å^3"

# 6) compute ion pairs from volume & salt molarity
echo "==> Step 6: Compute ion pairs for ${SALT_M} M"
N_PAIRS="$(ion_pairs_py "${V_A3}" "${SALT_M}")"; [[ -z "$N_PAIRS" ]] && N_PAIRS=0
echo "    Add ~${N_PAIRS} ${ION_POS} and ${N_PAIRS} ${ION_NEG}"

# 7) add salt ions and write final parm/rst
echo "==> Step 7: LEaP add salt"
cat > leap_add_salt.in <<EOF
set default PBRadii mbondi3
source ${FF_PROT}
source ${WATER_LEAP}
${UNIT}=loadpdb ${ROOT}_pqr.pdb
addIons ${UNIT} ${NEUT_ION} ${NEUT_COUNT}
${SOLVATE_CMD} ${UNIT} ${WATBOX} ${CUSHION_A}
addIonsRand ${UNIT} ${ION_POS} ${N_PAIRS}
addIonsRand ${UNIT} ${ION_NEG} ${N_PAIRS}
saveAmberParm ${UNIT} ${ROOT}_salt.prmtop ${ROOT}_salt.inpcrd
quit
EOF
run_cmd tleap -f leap_add_salt.in | tee leap_add_salt.log

# 8) cpptraj final PDB
echo "==> Step 8: Write final PDB with cpptraj"
cat > write_final_pdb.in <<EOF
parm ${ROOT}_salt.prmtop
trajin ${ROOT}_salt.inpcrd
autoimage
trajout ${ROOT}_final.pdb pdb
run
EOF
run_cmd cpptraj -i write_final_pdb.in | tee write_final_pdb.log

echo
echo "=== DONE ==="
[[ -f "$NEF" ]] && [[ -f "${ROOT}.RST" ]] && echo "  ${ROOT}.RST / ${ROOT}.DIP (from NEF)"
echo "  ${ROOT}_new.pdb (cleaned by pdb4amber)"
echo "  ${ROOT}_nocrys.pdb (CRYST1 removed)"
echo "  ${ROOT}_pqr.pdb / ${ROOT}.pqr"
echo "  ${ROOT}_neutral.prmtop/.inpcrd"
echo "  ${ROOT}_salt.prmtop/.inpcrd"
echo "  ${ROOT}_final.pdb"
"""
MINRELAX_SH = r"""#!/usr/bin/env bash
set -euo pipefail

# minrelax.sh — minimisation + heating + NPT relax
# Writes mdin files following the provided templates (1min/5md/2mdheat/3md/4md/6md/7md/8md/9md).
#
# Restraints logic:
#   - If --use_rst 1: write nmropt=1 + &wt/LISTOUT/ DISANG
#   - If --use_rst 0: write nmropt=0 and REMOVE &wt/LISTOUT/DISANG lines

unset DISPLAY || true
unset XAUTHORITY || true

WORKDIR=""
NP="8"
TEMP="300.0"
USE_RST="0"
RSTFILE=""
OUTDIR="02_minrelax"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --workdir) WORKDIR="$2"; shift 2;;
    --np) NP="$2"; shift 2;;
    --temp) TEMP="$2"; shift 2;;
    --use_rst) USE_RST="$2"; shift 2;;
    --rstfile) RSTFILE="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    *) echo "WARNING: unknown arg $1"; shift;;
  esac
done


[[ -z "$WORKDIR" ]] && { echo "ERROR: --workdir required"; exit 2; }

PREP="${WORKDIR}/01_preprocess"
PRMTOP="${PREP}/input_salt.prmtop"
INPCRD="${PREP}/input_salt.inpcrd"

[[ -f "$PRMTOP" ]] || { echo "ERROR: missing $PRMTOP. Run preprocessing first."; exit 2; }
[[ -f "$INPCRD" ]] || { echo "ERROR: missing $INPCRD. Run preprocessing first."; exit 2; }

# Convert key inputs to absolute paths so they remain valid after changing directory
PRMTOP="$(readlink -f "$PRMTOP")"
INPCRD="$(readlink -f "$INPCRD")"
PREP_ABS="$(readlink -f "$PREP")"

OUT="${WORKDIR}/${OUTDIR}"
mkdir -p "$OUT"
cd "$OUT"

# Prefer pmemd.MPI, but fall back to sander.MPI (AmberTools-only installs).
# Allow override from environment: AMBER_MPI_ENGINE / AMBER_SERIAL_ENGINE
PMEMD_MIN="${AMBER_SERIAL_ENGINE:-pmemd}"
PMEMD_MPI="${AMBER_MPI_ENGINE:-pmemd.MPI}"

# Engine selection by GPU request:
#   NGPUS=0  -> pmemd.MPI (CPU MPI)
#   NGPUS=1  -> pmemd.cuda (single-GPU)
#   NGPUS>=2 -> pmemd.cuda.MPI (multi-GPU MPI)
if [ "${NGPUS:-0}" -ge 2 ]; then
  PMEMD_MPI="pmemd.cuda.MPI"
elif [ "${NGPUS:-0}" -eq 1 ]; then
  PMEMD_MPI="pmemd.cuda"
else
  PMEMD_MPI="${AMBER_MPI_ENGINE:-pmemd.MPI}"
fi

# Auto-detect if the requested engine is missing
cmd_exists(){ command -v "$1" >/dev/null 2>&1; }

# Source AMBER if available on execute node (common on NMRbox)
if [ -f "/usr/software/amber/amber.sh" ]; then
  # shellcheck disable=SC1091
  source "/usr/software/amber/amber.sh"
fi

if ! cmd_exists "$PMEMD_MIN"; then
  if cmd_exists sander; then PMEMD_MIN="sander"; fi
fi

if ! cmd_exists "$PMEMD_MPI"; then
  if cmd_exists sander.MPI; then
    PMEMD_MPI="sander.MPI"
  elif cmd_exists pmemd; then
    # no MPI engine available; run serial pmemd for all stages
    PMEMD_MPI="pmemd"
  elif cmd_exists sander; then
    PMEMD_MPI="sander"
  fi
fi

cmd_exists "$PMEMD_MIN" || { echo "ERROR: No AMBER serial engine found (pmemd/sander). This execute node likely lacks AMBER. Use Condor requirements to target AMBER-ready nodes or ensure /usr/software/amber is available."; exit 1; }
cmd_exists "$PMEMD_MPI" || { echo "ERROR: No AMBER engine found (pmemd.MPI/sander.MPI/pmemd/sander)"; exit 1; }

if [ "${NP}" -lt 2 ]; then
  # Running with 1 proc: do NOT use MPI engine (pmemd.MPI requires >=2 ranks)
  if cmd_exists pmemd; then
    PMEMD_MPI="pmemd"
  elif cmd_exists sander; then
    PMEMD_MPI="sander"
  fi
fi


# Stage restraints file locally to avoid long-path issues in mdin namelists
DISANG_FILE=""
if [[ "$USE_RST" == "1" ]]; then
  if [[ -z "$RSTFILE" ]]; then
    RSTFILE="${PREP_ABS}/input.RST"
  fi
  [[ -f "$RSTFILE" ]] || { echo "ERROR: --use_rst=1 but RST file not found at $RSTFILE"; exit 2; }
  [[ -s "$RSTFILE" ]] || { echo "ERROR: RST file is empty: $RSTFILE"; exit 2; }

  DISANG_FILE="restraints.RST"
  cp -f "$RSTFILE" "$DISANG_FILE"
fi

write_nmr_tail() {
  # Append after &cntrl (i.e. after "/" line) ONLY when restraints are enabled
  if [[ "$USE_RST" == "1" ]]; then
    cat <<EOF
&wt TYPE="END", /
LISTOUT=POUT
DISANG=${DISANG_FILE}
EOF
  fi
}

write_nmropt_line() {
  if [[ "$USE_RST" == "1" ]]; then
    echo "  nmropt=1,"
  else
    echo "  nmropt=0,"
  fi
}



run_mpi() {
  local in="$1" out="$2" c="$3" r="$4" x="$5" ref="$6"
  # If selected engine is MPI and NP>1, launch with mpirun. Otherwise run directly.
  if [[ "$PMEMD_MPI" == *".MPI" ]] && command -v mpirun >/dev/null 2>&1 && [[ "${NP}" -gt 1 ]]; then
    mpirun -np "${NP}" "${PMEMD_MPI}" -O -i "$in" -o "$out" -p "$PRMTOP" -c "$c" -r "$r" -x "$x" -ref "$ref"
  else
    "${PMEMD_MPI}" -O -i "$in" -o "$out" -p "$PRMTOP" -c "$c" -r "$r" -x "$x" -ref "$ref"
  fi
}

# Reference for positional restraints
REFC="$INPCRD"

# ---------------- Stage 1: minimization of solvent (template 1min.in) ----------------
cat > 1min.in <<EOF
minimization of solvent
 &cntrl
  imin = 1, maxcyc = 1000,
  ncyc = 20, ntx = 1,
  ntwe = 0, ntwr = 500, ntpr = 50,
  ntc = 1, ntf = 1, ntb = 1, ntp = 0,
  cut = 10.0,
  ntr=1, restraintmask="!@H=", restraint_wt=100.0,
  ioutfm=1, ntxo=2,
$(write_nmropt_line)
 /
$(write_nmr_tail)
EOF

echo "==> Stage 1: 1min (serial pmemd)"
pmemd -O -i 1min.in -o 1min.out -p "$PRMTOP" -c "$INPCRD" -r 1min.rst7 -ref "$REFC"

# ---------------- Stage 3: heating (template 2mdheat.in) ----------------
cat > 2mdheat.in <<EOF
&cntrl
  imin = 0, nstlim = 100000, dt = 0.001,
  irest = 0, ntx = 1, ig = -1,
  tempi = 100.0, temp0 = ${TEMP},
  ntc = 2, ntf = 2, tol = 0.00001,
  ntwx = 10000, ntwe = 0, ntwr = 1000, ntpr = 1000,
  cut = 8.0, iwrap = 0,
  ntt = 3, gamma_ln = 1.0, ntb = 1, ntp = 0,
  nscm = 0,
  ntr=1, restraintmask="!@H=", restraint_wt=100.0,
$(write_nmropt_line)
  ioutfm=1, ntxo=2,
 /
EOF
# Only keep TEMP0 ramp when restraints are enabled (per user request: remove other lines otherwise)
if [[ "$USE_RST" == "1" ]]; then
  cat >> 2mdheat.in <<EOF
&wt TYPE="TEMP0", istep1=0, istep2=100000, value1=100.0, value2=${TEMP}, /
$(write_nmr_tail)
EOF
fi

echo "==> Stage 2: 2mdheat (pmemd.MPI; NP=${NP})"
run_mpi 2mdheat.in 2mdheat.out 1min.rst7 2mdheat.rst7 2mdheat.nc "$REFC"

# ---------------- Stage 4: NPT (template 3md.in) ----------------
cat > 3md.in <<EOF
 &cntrl
  imin = 0, nstlim = 100000, dt = 0.001,
  irest = 1, ntx = 5, ig = -1,
  temp0 = ${TEMP},
  ntc = 2, ntf = 2, tol = 0.00001,
  ntwx = 10000, ntwe = 0, ntwr = 1000, ntpr = 1000,
  cut = 8.0, iwrap = 0,
  ntt = 3, gamma_ln = 1.0, ntb = 2, ntp = 1, barostat = 2,
  nscm = 0,
  ntr=1, restraintmask="!@H=", restraint_wt=100.0,
  ioutfm=1, ntxo=2,
$(write_nmropt_line)
 /
$(write_nmr_tail)
EOF

echo "==> Stage 3: 3md (pmemd.MPI; NP=${NP})"
run_mpi 3md.in 3md.out 2mdheat.rst7 3md.rst7 3md.nc "$REFC"

# ---------------- Stage 5: NPT (template 4md.in) ----------------
cat > 4md.in <<EOF
 &cntrl
  imin = 0, nstlim = 100000, dt = 0.001,
  irest = 1, ntx = 5, ig = -1,
  temp0 = ${TEMP},
  ntc = 2, ntf = 2, tol = 0.00001,
  ntwx = 10000, ntwe = 0, ntwr = 1000, ntpr = 1000,
  cut = 8.0, iwrap = 0,
  ntt = 3, gamma_ln = 1.0, ntb = 2, ntp = 1, barostat = 2,
  nscm = 0,
  ntr=1, restraintmask="!@H=", restraint_wt=10.0,
  ioutfm=1, ntxo=2,
$(write_nmropt_line)
 /
$(write_nmr_tail)
EOF

echo "==> Stage 4: 4md (pmemd.MPI; NP=${NP})"
run_mpi 4md.in 4md.out 3md.rst7 4md.rst7 4md.nc "$REFC"

# ---------------- Stage 6: minimization excluding backbone (template 5min.in) ----------------
cat > 5min.in <<EOF
minimization of everything excluding backbone
 &cntrl
  imin = 1, maxcyc = 1000,
  ncyc = 30, ntx = 1,
  ntwe = 0, ntwr = 500, ntpr = 50,
  ntc = 1, ntf = 1, ntb = 1, ntp = 0,
  cut = 8.0,
  ntr=1, restraintmask="@CA,N,C", restraint_wt=10.0,
  ioutfm=1, ntxo=2,
$(write_nmropt_line)
 /
$(write_nmr_tail)
EOF

echo "==> Stage 5: 5min (serial pmemd)"
pmemd -O -i 5min.in -o 5min.out -p "$PRMTOP" -c 4md.rst7 -r 5min.rst7 -ref "$REFC"


# ---------------- Stage 6: NPT (template 6md.in) ----------------
cat > 6md.in <<EOF
&cntrl
  imin = 0, nstlim = 100000, dt = 0.001,
  irest = 0, ntx = 1, ig = -1,
  tempi = ${TEMP}, temp0 = ${TEMP},
  ntc = 2, ntf = 2, tol = 0.00001,
  ntwx = 10000, ntwe = 0, ntwr = 1000, ntpr = 1000,
  cut = 8.0, iwrap = 0,
  ntt = 3, gamma_ln = 1.0, ntb = 2, ntp = 1,
  nscm = 0, barostat = 2,
  ntr=1, restraintmask="@CA,N,C", restraint_wt=10.0,
  ioutfm=1, ntxo=2,
$(write_nmropt_line)
 /
$(write_nmr_tail)
EOF

echo "==> Stage 6: 6md (pmemd.MPI; NP=${NP})"
run_mpi 6md.in 6md.out 5min.rst7 6md.rst7 6md.nc "$REFC"

# ---------------- Stage 7: NPT (template 7md.in) ----------------
cat > 7md.in <<EOF
&cntrl
  imin = 0, nstlim = 100000, dt = 0.001,
  irest = 1, ntx = 5, ig = -1,
  temp0 = ${TEMP},
  ntc = 2, ntf = 2, tol = 0.00001,
  ntwx = 10000, ntwe = 0, ntwr = 1000, ntpr = 1000,
  cut = 8.0, iwrap = 0,
  ntt = 3, gamma_ln = 1.0, ntb = 2, ntp = 1,
  nscm = 0, barostat = 2,
  ntr=1, restraintmask="@CA,N,C", restraint_wt=1.0,
  ioutfm=1, ntxo=2,
$(write_nmropt_line)
 /
$(write_nmr_tail)
EOF

echo "==> Stage 7: 7md (pmemd.MPI; NP=${NP})"
run_mpi 7md.in 7md.out 6md.rst7 7md.rst7 7md.nc "$REFC"

# ---------------- Stage 8: NPT (template 8md.in) ----------------
cat > 8md.in <<EOF
&cntrl
  imin = 0, nstlim = 100000, dt = 0.001,
  irest = 1, ntx = 5, ig = -1,
  temp0 = ${TEMP},
  ntc = 2, ntf = 2, tol = 0.00001,
  ntwx = 10000, ntwe = 0, ntwr = 1000, ntpr = 1000,
  cut = 8.0, iwrap = 0,
  ntt = 3, gamma_ln = 1.0, ntb = 2, ntp = 1,
  nscm = 0, barostat = 2,
  ntr=1, restraintmask="@CA,N,C", restraint_wt=0.1,
  ioutfm=1, ntxo=2,
$(write_nmropt_line)
 /
$(write_nmr_tail)
EOF

echo "==> Stage 8: 8md (pmemd.MPI; NP=${NP})"
run_mpi 8md.in 8md.out 7md.rst7 8md.rst7 8md.nc "$REFC"

# ---------------- Stage 9: NPT (template 9md.in) ----------------
cat > 9md.in <<EOF
&cntrl
  imin = 0, nstlim = 100000, dt = 0.001,
  irest = 1, ntx = 5, ig = -1,
  temp0 = ${TEMP},
  ntc = 2, ntf = 2, tol = 0.00001,
  ntwx = 10000, ntwe = 0, ntwr = 1000, ntpr = 1000,
  cut = 8.0, iwrap = 0,
  ntt = 3, gamma_ln = 1.0, ntb = 2, ntp = 1,
  nscm = 1000, barostat = 2,
  ioutfm=1, ntxo=2,
$(write_nmropt_line)
 /
$(write_nmr_tail)
EOF

echo "==> Stage 9: 9md (pmemd.MPI; NP=${NP})"
run_mpi 9md.in 9md.out 8md.rst7 9md.rst7 9md.nc "$REFC"

echo "==> DONE. Final restart: ${OUT}/9md.rst7"
"""

PROD_SH = r"""#!/usr/bin/env bash
set -euo pipefail
# production.sh — Production MD (CPU or GPU) with optional NMR distance restraints (DISANG)
unset DISPLAY || true
unset XAUTHORITY || true

WORKDIR=""
NP="8"
NGPUS="0"
TEMP="300.0"
NSTLIM="250000"
DT="0.002"
USE_RST="0"
RSTFILE=""
OUTDIR="03_production"
START_RST=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --workdir) WORKDIR="$2"; shift 2;;
    --np) NP="$2"; shift 2;;
    --ngpus) NGPUS="$2"; shift 2;;
    --temp) TEMP="$2"; shift 2;;
    --nstlim) NSTLIM="$2"; shift 2;;
    --dt) DT="$2"; shift 2;;
    --use_rst) USE_RST="$2"; shift 2;;
    --rstfile) RSTFILE="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --start_rst) START_RST="$2"; shift 2;;
    *) echo "WARNING: unknown arg $1"; shift;;
  esac
done

[[ -z "$WORKDIR" ]] && { echo "ERROR: --workdir required"; exit 2; }

PREP="${WORKDIR}/01_preprocess"
PRMTOP="${PREP}/input_salt.prmtop"
[[ -f "$PRMTOP" ]] || { echo "ERROR: missing $PRMTOP. Run preprocessing first."; exit 2; }

if [[ -z "$START_RST" ]]; then
  START_RST="${WORKDIR}/02_minrelax/OUT/9md.rst7"
  [[ -f "$START_RST" ]] || START_RST="${WORKDIR}/02_minrelax/9md.rst7"
fi
[[ -f "$START_RST" ]] || { echo "ERROR: starting restart not found: $START_RST"; exit 2; }

PRMTOP="$(readlink -f "$PRMTOP")"
START_RST="$(readlink -f "$START_RST")"
PREP_ABS="$(readlink -f "$PREP")"

OUT="${WORKDIR}/${OUTDIR}"
mkdir -p "$OUT"
cd "$OUT"

cmd_exists(){ command -v "$1" >/dev/null 2>&1; }

if [ -f "/usr/software/amber/amber.sh" ]; then
  source "/usr/software/amber/amber.sh" || true
fi

# sanitize NP/NGPUS
if ! [[ "$NP" =~ ^[0-9]+$ ]] || [[ "$NP" -lt 1 ]]; then NP="1"; fi
if ! [[ "$NGPUS" =~ ^[0-9]+$ ]] || [[ "$NGPUS" -lt 0 ]]; then NGPUS="0"; fi
if [[ "$NGPUS" -ge 2 ]] && [[ "$NP" -lt 2 ]]; then NP="2"; fi

# choose engine
ENGINE=""
MODE="cpu"
if [[ "$NGPUS" -ge 2 ]]; then
  ENGINE="pmemd.cuda.MPI"; MODE="gpu_mpi"
elif [[ "$NGPUS" -eq 1 ]]; then
  ENGINE="pmemd.cuda"; MODE="gpu"
else
  ENGINE="pmemd.MPI"; MODE="cpu_mpi"
fi

# fallbacks
if ! cmd_exists "$ENGINE"; then
  if [[ "$ENGINE" == "pmemd.cuda.MPI" ]] && cmd_exists pmemd.cuda; then ENGINE="pmemd.cuda"; MODE="gpu"; fi
fi
if ! cmd_exists "$ENGINE"; then
  if [[ "$ENGINE" == "pmemd.MPI" ]] && cmd_exists sander.MPI; then ENGINE="sander.MPI"; fi
fi
if ! cmd_exists "$ENGINE"; then
  if [[ "$ENGINE" == "pmemd.cuda" ]] && cmd_exists pmemd; then ENGINE="pmemd"; MODE="cpu"; fi
fi
if ! cmd_exists "$ENGINE"; then
  if cmd_exists pmemd; then ENGINE="pmemd"; MODE="cpu"; fi
fi
cmd_exists "$ENGINE" || { echo "ERROR: No AMBER engine found on this node"; exit 1; }

# restraints
DISANG_FILE=""
if [[ "$USE_RST" == "1" ]]; then
  if [[ -z "$RSTFILE" ]]; then
    RSTFILE="${PREP_ABS}/input.RST"
  fi
  [[ -f "$RSTFILE" ]] || { echo "ERROR: --use_rst=1 but RST file not found at $RSTFILE"; exit 2; }
  [[ -s "$RSTFILE" ]] || { echo "ERROR: RST file is empty: $RSTFILE"; exit 2; }
  DISANG_FILE="restraints.RST"
  cp -f "$RSTFILE" "$DISANG_FILE"
fi

write_nmropt_line() { [[ "$USE_RST" == "1" ]] && echo "  nmropt=1," || echo "  nmropt=0,"; }

write_nmr_tail() {
  if [[ "$USE_RST" == "1" ]]; then
    cat <<EOF
&wt TYPE="END", /
LISTOUT=POUT
DISANG=${DISANG_FILE}
EOF
  fi
}

cat > prod.in <<EOF
Production MD (NPT)
 &cntrl
  imin=0, irest=1, ntx=5,
  nstlim=${NSTLIM}, dt=${DT},
  temp0=${TEMP},
  ntt=3, gamma_ln=1.0, ig=-1,
  ntc=2, ntf=2, tol=0.00001,
  ntwx=5000, ntpr=5000, ntwr=5000,
  ioutfm=1, ntxo=2,
  cut=8.0, iwrap=0,
  ntb=2, ntp=1, barostat=2, pres0=1.0, taup=1.0,
  nscm=1000,
$(write_nmropt_line)
 /
$(write_nmr_tail)
EOF

run_cmd() {
  local c="$1" o="$2" r="$3" x="$4"
  if [[ "$MODE" == "cpu_mpi" || "$MODE" == "gpu_mpi" ]]; then
    [[ "$NP" -lt 2 ]] && NP="2"
    mpirun -np "$NP" "$ENGINE" -O -i prod.in -o "$o" -p "$PRMTOP" -c "$c" -r "$r" -x "$x"
  else
    "$ENGINE" -O -i prod.in -o "$o" -p "$PRMTOP" -c "$c" -r "$r" -x "$x"
  fi
}

echo "==> Production: engine=$ENGINE mode=$MODE NP=$NP NGPUS=$NGPUS"
run_cmd "$START_RST" prod.out prod.rst7 prod.nc
echo "==> DONE. Final restart: ${OUT}/prod.rst7"
"""

ANALYSIS_CPPTRAJ_SH = r"""#!/usr/bin/env bash
set -euo pipefail

# analysis_cpptraj.sh
# Computes RMSD and RMSF using cpptraj.
#
# Flags:
#   --prmtop <file>   (required)
#   --traj <file>     (required)
#   --sel <cpptraj_mask>   (default ':1-999@C,CA,N,O' backbone heavy)
#   --ref <file>      optional reference (rst7/pdb); if not given uses first frame
#   --outdir <dir>    (default 04_analysis)

PRMTOP=""
TRAJ=""
SEL=":1-999@C,CA,N,O"
REF=""
OUTDIR="04_analysis"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --prmtop) PRMTOP="$2"; shift 2;;
    --traj) TRAJ="$2"; shift 2;;
    --sel) SEL="$2"; shift 2;;
    --ref) REF="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    *) echo "WARNING: unknown arg $1"; shift;;
  esac
done

[[ -f "$PRMTOP" ]] || { echo "ERROR: missing/invalid --prmtop"; exit 2; }
[[ -f "$TRAJ" ]]   || { echo "ERROR: missing/invalid --traj"; exit 2; }

mkdir -p "$OUTDIR"
cd "$OUTDIR"

REFLINE=""
if [[ -n "$REF" ]]; then
  [[ -f "$REF" ]] || { echo "ERROR: ref not found: $REF"; exit 2; }
  REFLINE="reference $REF"
fi

cat > cpptraj_analysis.in <<EOF
parm $PRMTOP
trajin $TRAJ
autoimage
$REFLINE

rms RMSD $SEL out rmsd.dat
atomicfluct RMSF $SEL byres out rmsf_byres.dat

run
quit
EOF

cpptraj -i cpptraj_analysis.in > cpptraj_analysis.log
echo "Wrote: $OUTDIR/rmsd.dat $OUTDIR/rmsf_byres.dat"
"""
PARSER_VIOLATIONS_PY = r"""#!/usr/bin/env python3
import argparse, re, csv
from pathlib import Path

def main():
    ap = argparse.ArgumentParser(description="Parse AMBER mdout for restraint violation listings (best with LISTOUT=POUT).")
    ap.add_argument("--mdout", required=True, help="AMBER output file (e.g. prod.out)")
    ap.add_argument("--out", default="violations.csv", help="Output CSV")
    args = ap.parse_args()

    mdout = Path(args.mdout)
    text = mdout.read_text(errors="replace").splitlines()

    step = None
    step_re = re.compile(r"\\bNSTEP\\s*=\\s*(\\d+)")
    keys = ("VIOL", "Violation", "DEVIAT", "deviation", "R1", "R2", "r =", "ENERGY")
    rows = []
    for line in text:
        m = step_re.search(line)
        if m:
            step = int(m.group(1))
            continue
        if any(k in line for k in keys):
            dev = None
            mdev = re.search(r"deviation\\s*=\\s*([-0-9.]+)", line, flags=re.IGNORECASE)
            if mdev:
                dev = float(mdev.group(1))
            rows.append({"step": step, "deviation": dev, "line": line.strip()})

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["step", "deviation", "line"])
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {out} with {len(rows)} rows")
    if len(rows) == 0:
        print("NOTE: For best results, run restrained MD with LISTOUT=POUT in mdin (DISANG).")

if __name__ == "__main__":
    main()
"""
PLOT_ANALYSIS_PY = r"""#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt

def read_two_cols(path):
    xs, ys = [], []
    for line in Path(path).read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#") or s.startswith("@"):
            continue
        parts = s.split()
        if len(parts) < 2:
            continue
        try:
            xs.append(float(parts[0]))
            ys.append(float(parts[1]))
        except Exception:
            continue
    return xs, ys

def read_rmsf_byres(path):
    # cpptraj atomicfluct byres typically: resnum rmsf
    return read_two_cols(path)

def plot_xy(x, y, xlabel, ylabel, title, outpng):
    if not x or not y:
        return False
    plt.figure()
    plt.plot(x, y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpng, dpi=200)
    plt.close()
    return True

def plot_violation_series(csv_path, outpng):
    import csv
    steps, devs = [], []
    with open(csv_path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            try:
                if row.get("deviation") in (None, "", "None"):
                    continue
                dev = float(row["deviation"])
                st = row.get("step")
                st = int(st) if st not in (None, "", "None") else None
                if st is None:
                    continue
                steps.append(st)
                devs.append(dev)
            except Exception:
                continue
    if not steps:
        return False
    plt.figure()
    plt.plot(steps, devs)
    plt.xlabel("Step")
    plt.ylabel("Deviation")
    plt.title("Restraint deviations (parsed)")
    plt.tight_layout()
    plt.savefig(outpng, dpi=200)
    plt.close()
    return True

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--analysis_dir", required=True, help="Directory containing rmsd.dat / rmsf_byres.dat / violations.csv")
    args = ap.parse_args()

    ad = Path(args.analysis_dir)
    rmsd = ad / "rmsd.dat"
    rmsf = ad / "rmsf_byres.dat"
    viol = ad / "violations.csv"

    if rmsd.exists():
        x, y = read_two_cols(rmsd)
        plot_xy(x, y, "Frame/Time", "RMSD (Å)", "RMSD", ad / "rmsd.png")
    if rmsf.exists():
        x, y = read_rmsf_byres(rmsf)
        plot_xy(x, y, "Residue", "RMSF (Å)", "RMSF by residue", ad / "rmsf_byres.png")
    if viol.exists():
        plot_violation_series(viol, ad / "violations.png")

if __name__ == "__main__":
    main()
"""


# ---------------- MEDOID UTILITIES (NMR ensemble PDB) ----------------


def _pdb_split_models(pdb_text):
    """Return list of model blocks (each block is list of ATOM/HETATM lines).
    If no MODEL records are present, returns a single model.
    """
    lines = pdb_text.splitlines(True)  # keep newlines
    models = []
    cur = []
    # in_model = False
    saw_model = False
    for ln in lines:
        rec = ln[0:6].strip()
        if rec == "MODEL":
            saw_model = True
            # in_model = True
            if cur:
                # flush any previous loose atoms (rare)
                models.append(cur)
                cur = []
            continue
        if rec == "ENDMDL":
            if cur:
                models.append(cur)
                cur = []
            # in_model = False
            continue
        if rec in ("ATOM", "HETATM"):
            cur.append(ln)
            continue
        # ignore other records for RMSD, but keep them implicitly by writing only ATOM/HETATM
    if cur:
        models.append(cur)
    if not saw_model:
        # treat whole file as single model: collect all ATOM/HETATM from original text
        models = [[ln for ln in lines if ln[0:6].strip() in ("ATOM", "HETATM")]]
    # drop empty models
    models = [m for m in models if m]
    return models


def _pdb_atom_key(line):
    """Create a stable atom identifier key across models."""
    # PDB columns: atom name 12-16, resname 17-20, chain 21, resid 22-26, insertion 26
    atom = line[12:16].strip()
    resn = line[17:20].strip()
    chain = line[21].strip()
    resid = line[22:26].strip()
    ins = line[26].strip()
    # alt = line[16].strip()  # altLoc
    # Ignore altLoc in key so we can intersect; alt handled by choosing first occurrence per key
    return (chain, resid, ins, resn, atom)


def _pdb_xyz(line):
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        return (x, y, z)
    except Exception:
        return None


def _is_protein_atom(line):
    # crude but effective: exclude waters/ions by resname; keep standard protein residues
    resn = line[17:20].strip().upper()
    if resn in ("HOH", "WAT", "SOL", "NA", "CL", "K", "CA", "MG", "ZN"):
        return False
    return True


def _choose_sel(atom_name, mode):
    atom_name = atom_name.upper()
    if mode == "ca":
        return atom_name == "CA"
    if mode == "bb":
        return atom_name in ("N", "CA", "C", "O")
    # heavy: exclude hydrogens
    return not atom_name.startswith("H")


def _kabsch_rmsd(P, Q):
    """RMSD after optimal superposition of P onto Q. P,Q shape (n,3)."""
    import numpy as np

    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    C = Pc.T @ Qc
    V, S, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    U = V @ D @ Wt
    P_rot = Pc @ U
    diff = P_rot - Qc
    return float(np.sqrt((diff * diff).sum() / P.shape[0]))


def _pick_medoid_model(pdb_path, sel_mode="bb", log_fn=None):
    """Return (medoid_index, models_lines). Index is 0-based."""
    txt = pathlib.Path(pdb_path).read_text()
    models = _pdb_split_models(txt)
    if len(models) <= 1:
        if log_fn:
            log_fn(
                "Medoid selection: input PDB contains a single model (no MODEL/ENDMDL ensemble detected); using as provided."
            )
        return 0, models

    try:
        import numpy as np  # noqa
    except Exception:
        if log_fn:
            log_fn(
                "WARNING: numpy not available; cannot compute medoid. Using model 1."
            )
        return 0, models

    # Build per-model coordinate dicts
    model_maps = []
    for mlines in models:
        d = {}
        for ln in mlines:
            if not _is_protein_atom(ln):
                continue
            atom = ln[12:16].strip()
            if not _choose_sel(atom, sel_mode):
                continue
            key = _pdb_atom_key(ln)
            if key in d:
                continue  # keep first occurrence
            xyz = _pdb_xyz(ln)
            if xyz is None:
                continue
            d[key] = xyz
        model_maps.append(d)

    common = set(model_maps[0].keys())
    for d in model_maps[1:]:
        common &= set(d.keys())

    if len(common) < 10:
        if log_fn:
            log_fn(
                f"WARNING: Only {len(common)} common atoms across models for medoid; using model 1."
            )
        return 0, models

    keys = sorted(common)
    coords = []
    for d in model_maps:
        coords.append([d[k] for k in keys])

    n = len(models)
    # Pairwise RMSD matrix (upper triangle)
    sums = [0.0] * n
    for i in range(n):
        Pi = coords[i]
        for j in range(i + 1, n):
            r = _kabsch_rmsd(Pi, coords[j])
            sums[i] += r
            sums[j] += r

    med = min(range(n), key=lambda i: sums[i])
    if log_fn:
        mean = sums[med] / (n - 1)
        log_fn(
            f"Medoid selection: chose model {med + 1}/{n} (mean pairwise RMSD {mean:.3f} Å) using '{sel_mode}' atoms."
        )
    return med, models


def _write_single_model_pdb(model_lines, out_path):
    with open(out_path, "w") as f:
        for ln in model_lines:
            if ln[0:6].strip() in ("ATOM", "HETATM"):
                f.write(ln)
        f.write("END\n")


class MdWorkflowDialog(CcpnDialog):
    def __init__(self, parent=None, mainWindow=None, **kwds):
        super().__init__(
            parent=parent, setLayout=True, windowTitle="NMR–MD Workflow (AMBER)", **kwds
        )

        self.proc = None

        self.mainFrame = Frame(self, setLayout=True, grid=(0, 0))
        self.tabs = QTabWidget()
        self.mainFrame.getLayout().addWidget(self.tabs)

        self._createInputTab()
        self._createPreprocessTab()
        self._createMinRelaxTab()
        self._createProductionTab()
        self._createAnalysisTab()

        self.logFrame = Frame(self, setLayout=True, grid=(1, 0))
        Label(self.logFrame, text="Command log:", grid=(0, 0))
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.logFrame.getLayout().addWidget(self.log)

        self._log("GUI initialised")

    # ---------------- utilities ----------------

    def _log(self, msg):
        ts = datetime.datetime.now().strftime("%H:%M:%S")
        self.log.appendPlainText(f"[{ts}] {msg}")

    def _getWorkDir(self) -> str:
        """Return the workflow work directory from the Input tab."""
        try:
            wd = (self.workDir.get() or "").strip()
        except Exception:
            # Fallback for Qt widgets
            try:
                wd = (self.workDir.text() or "").strip()
            except Exception:
                wd = ""
        return wd

    def _run(self, name, cmd, cwd):
        if self.proc and self.proc.state() != QProcess.NotRunning:
            showWarning("Running", "Another process is already running")
            return

        self.proc = QProcess(self)
        self.proc.setWorkingDirectory(cwd)
        self.proc.setProcessChannelMode(QProcess.MergedChannels)
        self.proc.readyRead.connect(
            lambda: self._log(bytes(self.proc.readAll()).decode(errors="replace"))
        )
        self.proc.finished.connect(lambda c, s: self._log(f"[DONE:{name}]"))

        self._log(f"[RUN:{name}]")
        self._log(cmd)
        self.proc.start(SHELL, ["-lc", cmd])

    def _runCommand(self, cmd, tag="CMD", cwd=None):
        """Run a command (list or string) via the same QProcess logger used by _run()."""
        if isinstance(cmd, (list, tuple)):
            parts = []
            for a in cmd:
                a = str(a)
                # allow simple $VARS to expand in bash -lc
                if re.fullmatch(r"\$[A-Za-z_][A-Za-z0-9_]*", a):
                    parts.append(a)
                else:
                    parts.append(shlex.quote(a))
            cmdStr = " ".join(parts)
        else:
            cmdStr = str(cmd)

        self._run(tag, cmdStr, cwd=cwd)

    def _ensure_script(self, wd, name, content):
        scripts_dir = os.path.join(wd, "scripts")
        os.makedirs(scripts_dir, exist_ok=True)
        script_path = os.path.join(scripts_dir, name)
        with open(script_path, "w") as fh:
            fh.write(content)
        os.chmod(script_path, stat.S_IRWXU)
        return script_path

    def _load_input_json(self, wd):
        p = os.path.join(wd, "input.json")
        if os.path.exists(p):
            try:
                with open(p) as fh:
                    return json.load(fh)
            except Exception:
                return None

    def _shape_flag_from_ui(self):
        """Map UI box-shape label to script flag."""
        txt = (self._txt(self.boxShape) or "").strip().lower()
        if txt.startswith("truncated") or txt == "oct":
            return "oct"
        if txt.startswith("cube"):
            return "cube"
        # Fallback: allow future shapes to pass through
        return txt or "oct"

    # ---------------- INPUT TAB ----------------
    def _idx(self, widget, default=0):
        """Return a robust 0-based selected index for CCPN widgets.

        CCPN widgets differ slightly across versions; PulldownList/RadioButtons
        may expose different APIs depending on CCPN build.
        """
        if widget is None:
            return default

        # Common CCPN widget APIs
        for name in (
            "getIndex",  # some CCPN widgets
            "getSelectedIndex",  # some RadioButtons variants
            "currentIndex",  # Qt-style
        ):
            if hasattr(widget, name) and callable(getattr(widget, name)):
                try:
                    v = getattr(widget, name)()
                    if v is None:
                        continue
                    return int(v)
                except Exception:
                    pass

        # Some widgets store index as an attribute
        for name in ("index", "_index"):
            if hasattr(widget, name):
                try:
                    return int(getattr(widget, name))
                except Exception:
                    pass

        # Last resort: infer from text, if possible
        try:
            txt = self._txt(widget, default=None)
            if txt is None:
                return default
            if hasattr(widget, "texts") and isinstance(widget.texts, (list, tuple)):
                if txt in widget.texts:
                    return int(widget.texts.index(txt))
        except Exception:
            pass

        return default

    def _txt(self, widget, default=""):
        """Return selected text from CCPN widgets robustly."""
        if widget is None:
            return default
        for name in ("getText", "get", "currentText", "text"):
            if hasattr(widget, name) and callable(getattr(widget, name)):
                try:
                    v = getattr(widget, name)()
                    if v is None:
                        continue
                    return str(v)
                except Exception:
                    pass
        # Try mapping index to texts list
        try:
            i = self._idx(widget, default=None)
            if isinstance(i, int):
                texts = None
                for tname in ("texts", "_texts"):
                    if hasattr(widget, tname):
                        texts = getattr(widget, tname)
                        break
                if texts and 0 <= i < len(texts):
                    return str(texts[i])
        except Exception:
            pass
        return default
        # Most CCPN widgets (RadioButtons, CheckBox with texts, etc.)
        for name in ("getIndex", "getSelectedIndex", "currentIndex", "getCurrentIndex"):
            if hasattr(widget, name) and callable(getattr(widget, name)):
                try:
                    v = getattr(widget, name)()
                    if isinstance(v, int):
                        return v
                except Exception:
                    pass
        # PulldownList / ComboBox-like widgets often have get() or currentText()
        for name in ("get", "currentText", "text"):
            if hasattr(widget, name) and callable(getattr(widget, name)):
                try:
                    v = getattr(widget, name)()
                    # Some CCPN widgets return the selected text
                    if isinstance(v, int):
                        return v
                    if isinstance(v, str):
                        # Map text back to an index if we can
                        texts = None
                        for tname in ("texts", "_texts"):
                            if hasattr(widget, tname):
                                texts = getattr(widget, tname)
                                break
                        if (
                            texts is None
                            and hasattr(widget, "getTexts")
                            and callable(getattr(widget, "getTexts"))
                        ):
                            try:
                                texts = widget.getTexts()
                            except Exception:
                                texts = None
                        if texts:
                            try:
                                return list(texts).index(v)
                            except ValueError:
                                return default
                except Exception:
                    pass
        return default

    def _createInputTab(self):
        tab = Frame(setLayout=True)
        self.tabs.addTab(tab, "Input")

        row = 0
        Label(tab, "Working directory:", grid=(row, 0))
        self.workDir = Entry(tab, grid=(row, 1))
        Button(tab, "Browse", grid=(row, 2), callback=self._browseDir)
        row += 1

        Label(tab, "Input PDB file:", grid=(row, 0))
        self.pdbFile = Entry(tab, grid=(row, 1))
        Button(tab, "Browse", grid=(row, 2), callback=self._browsePdb)
        row += 1

        Label(tab, "Input NEF file (optional):", grid=(row, 0))
        self.nefFile = Entry(tab, grid=(row, 1))
        Button(tab, "Browse", grid=(row, 2), callback=self._browseNef)
        row += 1

        Label(tab, "pH:", grid=(row, 0))
        self.ph = FloatEntry(tab, text=7.0, grid=(row, 1))
        row += 1

        Label(tab, "Temperature (K):", grid=(row, 0))
        self.temp = FloatEntry(tab, text=300.0, grid=(row, 1))
        row += 1

        Label(tab, "Salt concentration (M):", grid=(row, 0))
        self.salt = FloatEntry(tab, text=0.150, grid=(row, 1))
        row += 1

        Button(tab, "Save input", grid=(row, 1), callback=self._saveInput)

    def _browseDir(self):
        d = QFileDialog.getExistingDirectory(self, "Select working directory")
        if d:
            self.workDir.setText(d)

    def _browsePdb(self):
        f, _ = QFileDialog.getOpenFileName(self, "Select PDB", filter="PDB (*.pdb)")
        if f:
            self.pdbFile.setText(f)

    def _browseNef(self):
        f, _ = QFileDialog.getOpenFileName(self, "Select NEF", filter="NEF (*.nef)")
        if f:
            self.nefFile.setText(f)

    def _saveInput(self):
        wd = self.workDir.get()
        if not wd or not self.pdbFile.get():
            showWarning("Missing input", "Working directory and PDB file are required")
            return

        os.makedirs(os.path.join(wd, "source"), exist_ok=True)
        shutil.copy2(self.pdbFile.get(), os.path.join(wd, "source", "input.pdb"))

        # If requested, choose medoid model from an NMR ensemble PDB and overwrite staged input.pdb
        try:
            if hasattr(self, "ensembleMode") and self._idx(self.ensembleMode) == 1:
                staged = os.path.join(wd, "source", "input.pdb")
                # Keep a copy of the ensemble
                ensemble_copy = os.path.join(wd, "source", "input_ensemble.pdb")
                if not os.path.exists(ensemble_copy):
                    shutil.copy2(staged, ensemble_copy)
                # Selection mode
                sel = "bb"
                if hasattr(self, "medoidAtoms"):
                    idx = self._idx(self.medoidAtoms)
                    sel = "bb" if idx == 0 else ("ca" if idx == 1 else "heavy")
                med_i, models = _pick_medoid_model(
                    staged, sel_mode=sel, log_fn=self._log
                )
                med_path = os.path.join(wd, "source", "input_medoid.pdb")
                _write_single_model_pdb(models[med_i], med_path)
                shutil.copy2(med_path, staged)
        except Exception as e:
            self._log(
                f"WARNING: Medoid selection failed; using provided PDB. Details: {e}"
            )

        nef = self.nefFile.get()
        if nef:
            shutil.copy2(nef, os.path.join(wd, "source", "input.nef"))

        data = {
            "pdb": "source/input.pdb",
            "nef": "source/input.nef" if nef else None,
            "ph": float(self.ph.get()),
            "temperature": float(self.temp.get()),
            "salt_M": float(self.salt.get()),
        }

        with open(os.path.join(wd, "input.json"), "w") as fh:
            json.dump(data, fh, indent=2)

        self._log("Input saved")

    # ---------------- PREPROCESS TAB ----------------

    def _createPreprocessTab(self):
        tab = Frame(setLayout=True)
        self.tabs.addTab(tab, "Pre-processing")

        row = 0
        Label(tab, "Protein force field:", grid=(row, 0))
        self.ffProt = PulldownList(tab, texts=["ff19SB", "ff14SB"], grid=(row, 1))
        row += 1

        Label(tab, "Water model:", grid=(row, 0))
        self.ffWat = PulldownList(tab, texts=["tip3p", "spce", "opc"], grid=(row, 1))
        row += 1

        Label(tab, "Protonation method:", grid=(row, 0))
        self.protMethod = RadioButtons(
            tab,
            texts=["pdb2pqr + PropKa", "pypka (inactive)"],
            selectedInd=0,
            grid=(row, 1),
        )
        row += 1

        # Ensemble handling (NMR multi-model PDB)
        Label(tab, "Ensemble PDB:", grid=(row, 0))
        self.ensembleMode = RadioButtons(
            tab,
            texts=["Use first model", "Choose medoid (best representative)"],
            selectedInd=0,
            grid=(row, 1),
        )
        row += 1

        Label(tab, "Medoid RMSD atoms:", grid=(row, 0))
        self.medoidAtoms = PulldownList(
            tab,
            texts=["Backbone (C,CA,N,O)", "C-alpha (CA)", "Protein heavy atoms"],
            grid=(row, 1),
        )
        self.medoidAtoms.setIndex(0)
        row += 1

        Label(tab, "Box shape:", grid=(row, 0))
        self.boxShape = PulldownList(
            tab, texts=["Truncated Octahedron", "Cube"], grid=(row, 1)
        )
        row += 1

        Label(tab, "Box buffer (Å):", grid=(row, 0))
        self.box = FloatEntry(tab, text=12.0, grid=(row, 1))
        row += 1

        Label(tab, "nmrmd venv (path):", grid=(row, 0))
        self.venvPath = Entry(tab, text=DEFAULT_VENV, grid=(row, 1))
        row += 1

        Button(tab, "Run pre-processing", grid=(row, 1), callback=self._runPreprocess)

    def _runPreprocess(self):
        wd = self.workDir.get()
        if not wd:
            showWarning("Missing input", "Select working directory first")
            return

        # Ensure inputs staged (if user didn't click Save input)
        if self.pdbFile.get():
            os.makedirs(os.path.join(wd, "source"), exist_ok=True)
            shutil.copy2(self.pdbFile.get(), os.path.join(wd, "source", "input.pdb"))

        # If requested, choose medoid model from an NMR ensemble PDB and overwrite staged input.pdb
        try:
            if hasattr(self, "ensembleMode") and self._idx(self.ensembleMode) == 1:
                staged = os.path.join(wd, "source", "input.pdb")
                # Keep a copy of the ensemble
                ensemble_copy = os.path.join(wd, "source", "input_ensemble.pdb")
                if not os.path.exists(ensemble_copy):
                    shutil.copy2(staged, ensemble_copy)
                # Selection mode
                sel = "bb"
                if hasattr(self, "medoidAtoms"):
                    idx = self._idx(self.medoidAtoms)
                    sel = "bb" if idx == 0 else ("ca" if idx == 1 else "heavy")
                med_i, models = _pick_medoid_model(
                    staged, sel_mode=sel, log_fn=self._log
                )
                med_path = os.path.join(wd, "source", "input_medoid.pdb")
                _write_single_model_pdb(models[med_i], med_path)
                shutil.copy2(med_path, staged)
        except Exception as e:
            self._log(
                f"WARNING: Medoid selection failed; using provided PDB. Details: {e}"
            )
        nef = self.nefFile.get()
        if nef:
            os.makedirs(os.path.join(wd, "source"), exist_ok=True)
            shutil.copy2(nef, os.path.join(wd, "source", "input.nef"))

        prep_script = self._ensure_script(wd, "prep_pdb.sh", PREP_PDB_SH)

        venv = (self.venvPath.get() or "").strip() or DEFAULT_VENV

        if self._idx(self.protMethod) != 0:
            self._log("WARNING: pypka is inactive; using pdb2pqr + PropKa")

        cmd = (
            "set -e; "
            f"source '{venv}/bin/activate' && "
            f"source '{AMBER_SH}' && "
            f"'{prep_script}' "
            f"--workdir '{wd}' "
            f"--ph {float(self.ph.get())} "
            f"--salt {float(self.salt.get())} "
            f"--protein_ff {self._txt(self.ffProt)} "
            f"--water {self._txt(self.ffWat)} "
            f"--box {float(self.box.get())} "
            f"--shape {self._shape_flag_from_ui()}"
        )

        self._run("PREPROCESS", cmd, wd)

    # ---------------- MIN/RELAX TAB ----------------

    def _createMinRelaxTab(self):
        tab = Frame(setLayout=True)
        self.tabs.addTab(tab, "Minimisation/Relaxation")

        row = 0
        Label(tab, "Include NMR restraints:", grid=(row, 0))
        self.useRst = RadioButtons(
            tab, texts=["No", "Yes"], selectedInd=0, grid=(row, 1)
        )
        row += 1

        Label(tab, "CPUs for relaxation (pmemd.MPI):", grid=(row, 0))
        self.nCpu = Entry(tab, text="8", grid=(row, 1))
        row += 1

        Label(tab, "Output folder name:", grid=(row, 0))
        self.minrelaxOut = Entry(tab, text="02_minrelax", grid=(row, 1))
        row += 1

        Button(
            tab,
            "Run minimisation + relaxation",
            grid=(row, 1),
            callback=self._runMinRelax,
        )
        row += 1

        # ---------------- HTCondor (HPC) for Min/Relax ----------------
        Label(tab, "HTCondor (HPC) — Min/Relax:", grid=(row, 0))
        row += 1

        Label(tab, "HTCondor CPUs:", grid=(row, 0))
        self.condorMinRelaxCpus = Entry(tab, text="16", grid=(row, 1))
        row += 1

        Label(tab, "HTCondor GPUs:", grid=(row, 0))
        self.condorMinRelaxGpus = Entry(tab, text="0", grid=(row, 1))
        row += 1

        Label(tab, "HTCondor memory:", grid=(row, 0))
        self.condorMinRelaxMem = Entry(tab, text="32 GB", grid=(row, 1))
        row += 1

        Button(
            tab,
            "Generate HTCondor scripts",
            grid=(row, 0),
            callback=self._condorMinRelaxGenerate,
        )
        Button(
            tab,
            "Submit to HTCondor",
            grid=(row, 1),
            callback=self._condorMinRelaxSubmit,
        )
        Button(
            tab, "Check job status", grid=(row, 2), callback=self._condorMinRelaxStatus
        )

    def _runMinRelax(self):
        wd = self.workDir.get()
        if not wd:
            showWarning("Missing input", "Select working directory first")
            return

        # Prefer temperature stored from input tab
        j = self._load_input_json(wd)
        temp = float(self.temp.get())
        if j and "temperature" in j and j["temperature"] is not None:
            try:
                temp = float(j["temperature"])
            except Exception:
                pass

        min_script = self._ensure_script(wd, "minrelax.sh", MINRELAX_SH)

        try:
            np = int(str(self.nCpu.get()).strip())
            if np < 1:
                np = 1
        except Exception:
            np = 1
            self._log("WARNING: invalid CPU count; using 1")

        use_rst = 1 if self._idx(self.useRst) == 1 else 0
        outdir = (self.minrelaxOut.get() or "02_minrelax").strip()

        cmd = (
            "set -e; "
            f"source '{AMBER_SH}' && "
            f"'{min_script}' "
            f"--workdir '{wd}' --np {np} --temp {temp} --use_rst {use_rst} --outdir '{outdir}'"
        )

        self._run("MINRELAX", cmd, wd)

    # ---------------- HTCondor (Min/Relax) ----------------

    def _condor_minrelax_dir(self, wd):
        return os.path.join(wd, "02_minrelax", "htcondor_minrelax")

    def _condorMinRelaxGenerate(self):
        maxrt = "259200"
        wd = self.workDir.get()
        if not wd:
            showWarning("Missing input", "Select working directory first")
            return

        condor_dir = self._condor_minrelax_dir(wd)
        os.makedirs(condor_dir, exist_ok=True)
        os.makedirs(os.path.join(condor_dir, "01_preprocess"), exist_ok=True)
        os.makedirs(os.path.join(condor_dir, "logs"), exist_ok=True)

        # Find topology + starting coordinates from preprocessing
        prep_dir = os.path.join(wd, "01_preprocess")

        prmtop = os.path.join(prep_dir, "input_salt.prmtop")
        if not os.path.isfile(prmtop):
            # fallback: any prmtop
            cands = [
                os.path.join(prep_dir, f)
                for f in os.listdir(prep_dir)
                if f.endswith(".prmtop")
            ]
            prmtop = cands[0] if cands else None
        if not prmtop or not os.path.isfile(prmtop):
            showWarning("Missing files", f"Cannot find prmtop in {prep_dir}")
            return

        # coordinate candidates
        coord = None
        for fn in (
            "input_salt.rst7",
            "input_salt.inpcrd",
            "input_salt.crd",
            "input_salt.rst",
            "input_salt.rst7.gz",
            "input_salt.inpcrd.gz",
            "input_salt.crd.gz",
            "input_salt.rst.gz",
        ):
            p = os.path.join(prep_dir, fn)
            if os.path.isfile(p):
                coord = p
                break
        if coord is None:
            # fallback: any rst7/inpcrd/crd
            for f in os.listdir(prep_dir):
                if f.endswith(
                    (
                        ".rst7",
                        ".inpcrd",
                        ".crd",
                        ".rst",
                        ".rst7.gz",
                        ".inpcrd.gz",
                        ".crd.gz",
                        ".rst.gz",
                    )
                ):
                    coord = os.path.join(prep_dir, f)
                    break
        if coord is None:
            showWarning(
                "Missing files", f"Cannot find starting coordinates in {prep_dir}"
            )
            return

        use_rst = 1 if self._idx(self.useRst) == 1 else 0
        rstfile = os.path.join(prep_dir, "input.RST")
        if use_rst and not os.path.isfile(rstfile):
            # fallback name
            alt = os.path.join(prep_dir, "input.RST.txt")
            if os.path.isfile(alt):
                rstfile = alt
            else:
                showWarning(
                    "Missing files",
                    f"Restraints requested but cannot find {prep_dir}/input.RST",
                )
                return

        # Copy inputs into the condor package preserving expected structure
        shutil.copy2(
            prmtop, os.path.join(condor_dir, "01_preprocess", os.path.basename(prmtop))
        )
        shutil.copy2(
            coord, os.path.join(condor_dir, "01_preprocess", os.path.basename(coord))
        )
        if use_rst:
            shutil.copy2(
                rstfile,
                os.path.join(condor_dir, "01_preprocess", os.path.basename(rstfile)),
            )

        # Write scripts
        minrelax_path = os.path.join(condor_dir, "minrelax.sh")
        pathlib.Path(minrelax_path).write_text(MINRELAX_SH)
        os.chmod(minrelax_path, 0o755)

        runner = os.path.join(condor_dir, "run_pipeline.sh")
        runner_text = r"""#!/bin/bash
set -e

# Determine NPROCS deterministically (avoid nproc/thread-count oversubscription)
ARG_NP="${1:-}"
if [[ "$ARG_NP" =~ ^[0-9]+$ ]] && [ "$ARG_NP" -gt 0 ]; then
  NPROCS="$ARG_NP"
elif [ -n "${NPROCS:-}" ] && [[ "${NPROCS}" =~ ^[0-9]+$ ]] && [ "${NPROCS}" -gt 0 ]; then
  NPROCS="${NPROCS}"
elif [ -n "${_CONDOR_SLOT_CPUS:-}" ] && [[ "${_CONDOR_SLOT_CPUS}" =~ ^[0-9]+$ ]] && [ "${_CONDOR_SLOT_CPUS}" -gt 0 ]; then
  NPROCS="${_CONDOR_SLOT_CPUS}"
else
  NPROCS="16"
fi
# Hard safety: never run pmemd.MPI with <2 ranks; force 16 as requested
if [ "$NPROCS" -lt 2 ]; then
  NPROCS="16"
fi
export NPROCS
# ---- helpers ----
cmd_exists(){ command -v "$1" >/dev/null 2>&1; }

OUT="OUT"
mkdir -p "$OUT"
echo "Job started: $(date) on $(hostname)" > "$OUT/_job_meta.txt"

# Try to make Amber available (works on NMRbox nodes)
if [ -r /usr/software/amber/amber.sh ]; then
  source /usr/software/amber/amber.sh || true
fi

ENGINE="pmemd.MPI"
# Engine selection by GPU request:
#   NGPUS=0  -> pmemd.MPI (CPU MPI)
#   NGPUS=1  -> pmemd.cuda (single-GPU)
#   NGPUS>=2 -> pmemd.cuda.MPI (multi-GPU MPI)
if [ "${NGPUS:-0}" -ge 2 ]; then
  ENGINE="pmemd.cuda.MPI"
elif [ "${NGPUS:-0}" -eq 1 ]; then
  ENGINE="pmemd.cuda"
else
  ENGINE="pmemd.MPI"
fi
LAUNCH="mpirun"

# If missing, try environment modules + common roots (NMRbox pool differences)
if ! cmd_exists "$ENGINE" || ! cmd_exists "$LAUNCH"; then
  if [ -r /etc/profile.d/modules.sh ]; then
    . /etc/profile.d/modules.sh || true
    for m in amber/24 amber/23 amber/22 amber/20 amber openmpi openmpi/4 openmpi/3; do
      module load "$m" >/dev/null 2>&1 || true
    done
  fi
  for base in /usr/software /opt /usr/local /usr; do
    for dir in "$base"/amber* "$base"/Amber*; do
      [ -d "$dir/bin" ] && PATH="$dir/bin:$PATH"
    done
  done
fi

echo "PATH=$PATH" >> "$OUT/_job_meta.txt"

# Pick best available AMBER engine on this execute node
if cmd_exists pmemd.MPI; then
  export AMBER_MPI_ENGINE="pmemd.MPI"
elif cmd_exists sander.MPI; then
  export AMBER_MPI_ENGINE="sander.MPI"
elif cmd_exists pmemd; then
  export AMBER_MPI_ENGINE="pmemd"
elif cmd_exists sander; then
  export AMBER_MPI_ENGINE="sander"
fi
echo "AMBER_MPI_ENGINE=${AMBER_MPI_ENGINE:-unset}" >> "$OUT/_job_meta.txt"
echo "AMBERHOME=${AMBERHOME:-unset}" >> "$OUT/_job_meta.txt"
for b in pmemd pmemd.MPI sander sander.MPI; do
  if cmd_exists "$b"; then
    echo "FOUND:$b=$(command -v $b)" >> "$OUT/_job_meta.txt"
  else
    echo "MISSING:$b" >> "$OUT/_job_meta.txt"
  fi
done


if ! cmd_exists "$ENGINE"; then
  echo "FATAL: $ENGINE not found on node." | tee -a "$OUT/_job_meta.txt"
  exit 127
fi
# Only require mpirun when using an MPI engine
if [[ "$ENGINE" == *".MPI" ]]; then
  if ! cmd_exists "$LAUNCH"; then
    echo "FATAL: mpirun not found on node." | tee -a "$OUT/_job_meta.txt"
    exit 127
  fi
fi

echo "ENGINE_PREF=$ENGINE" | tee -a "$OUT/_job_meta.txt"
if [[ "$ENGINE" == *".MPI" ]]; then
  echo "MPIRUN=$LAUNCH -np $NPROCS" >> "$OUT/_job_meta.txt"
else
  echo "MPIRUN=none" >> "$OUT/_job_meta.txt"
fi
USE_RST="${USE_RST:-0}"

# Run the same minrelax.sh protocol but inside this condor sandbox.
# Inputs are staged under ./01_preprocess/
mkdir -p 01_preprocess
# If Condor transferred files without subdirectories, stage them where minrelax.sh expects
for f in input_salt.prmtop input_salt.inpcrd input.RST; do
  if [ -f "./$f" ] && [ ! -f "./01_preprocess/$f" ]; then
    cp -f "./$f" "./01_preprocess/$f"
  fi
done

bash ./minrelax.sh --workdir "." --np "$NPROCS" --temp "${TEMP:-300.0}" --use_rst "$USE_RST" --outdir "OUT"

echo "Job finished: $(date)" >> "$OUT/_job_meta.txt"
"""
        pathlib.Path(runner).write_text(runner_text)
        os.chmod(runner, 0o755)

        # Submit file
        try:
            cpus = int(str(self.condorMinRelaxCpus.get()).strip())
            gpus = (
                int(getattr(self, "condorMinRelaxGpus", None).get() or 0)
                if hasattr(self, "condorMinRelaxGpus")
                else 0
            )
            if cpus < 1:
                cpus = 1
        except Exception:
            cpus = 1
        mem = (self.condorMinRelaxMem.get() or "32 GB").strip()
        submit_path = os.path.join(condor_dir, "amber_pmemd_minrelax.submit")
        submit_lines = [
            "universe                = vanilla",
            "initialdir             = .",
            "executable              = /bin/bash",
            "arguments               = run_pipeline.sh",
            "",
            "# CPU resources only",
            f"request_cpus            = {cpus}",
            f"request_gpus            = {gpus}",
            f"request_memory          = {mem}",
            f"+MaxRuntime             = {maxrt}",
            "",
            "# Pass NPROCS to wrapper (matches request_cpus)",
            'environment             = "NPROCS=$(request_cpus) NGPUS=$(request_gpus)"',
            "getenv                  = True",
            "",
            "transfer_executable     = false",
            "should_transfer_files   = YES",
            "when_to_transfer_output = ON_EXIT",
            "",
            "# Inputs",
            "transfer_input_files    = run_pipeline.sh, minrelax.sh, 01_preprocess/"
            + os.path.basename(prmtop)
            + ", 01_preprocess/"
            + os.path.basename(coord)
            + (", 01_preprocess/" + os.path.basename(rstfile) if use_rst else ""),
            "transfer_output_files   = OUT",
            "",
            "# Optional: only run on Amber-ready nodes",
            "# requirements            = (HasAMBER == true)",
            "",
            'requirements            = (AMBER =!= undefined) && (AMBERTOOLS =!= undefined) && regexp("nmrbox.org", Machine)',
            "log                     = condor.log",
            "output                  = condor.out",
            "error                   = condor.err",
            "",
            "queue 1",
            "",
        ]
        pathlib.Path(submit_path).write_text("\n".join(submit_lines))

        self._log("[CONDOR:MINRELAX] Generated HTCondor scripts")
        self._log(f"[CONDOR:MINRELAX]  dir   = {condor_dir}")
        self._log(f"[CONDOR:MINRELAX]  submit= {os.path.basename(submit_path)}")
        self._log(f"[CONDOR:MINRELAX]  exe   = {os.path.basename(runner)}")
        self._log(
            f"[CONDOR:MINRELAX]  inputs= {os.path.basename(prmtop)}, {os.path.basename(coord)}"
            + (f", {os.path.basename(rstfile)}" if use_rst else "")
        )

    def _condorMinRelaxSubmit(self):
        wd = self.workDir.get()
        if not wd:
            showWarning("Missing input", "Select working directory first")
            return
        condor_dir = self._condor_minrelax_dir(wd)
        submit_path = os.path.join(condor_dir, "amber_pmemd_minrelax.submit")
        if not os.path.isfile(submit_path):
            showWarning("Missing files", "Generate HTCondor scripts first")
            return

        cmd = f"condor_submit '{submit_path}'"
        self._run("CONDOR_MINRELAX", cmd, condor_dir)

        # Best-effort parse cluster id from condor_submit stdout by reading condor.log later
        self.lastCondorMinRelaxSubmit = datetime.datetime.now().isoformat()

    def _condorMinRelaxStatus(self):
        wd = self.workDir.get()
        if not wd:
            showWarning("Missing input", "Select working directory first")
            return
        condor_dir = self._condor_minrelax_dir(wd)

        # Try to read last cluster from submit output stored in condor.log
        # If not found, still show condor_q for the user.
        condor_log = os.path.join(condor_dir, "condor.log")
        cluster = None
        if os.path.isfile(condor_log):
            try:
                t = pathlib.Path(condor_log).read_text(errors="replace")
                m = re.search(r"submitted to cluster\s+(\d+)", t)
                if m:
                    cluster = m.group(1)
            except Exception:
                pass

        if cluster:
            self._run("CONDOR_Q", f"condor_q {cluster}", condor_dir)
            self._run("CONDOR_HIST", f"condor_history -limit 5 {cluster}", condor_dir)
        else:
            self._run("CONDOR_Q", "condor_q -submitter $USER", condor_dir)
            self._log(
                "[CONDOR:MINRELAX] Tip: open condor.log to find the cluster id if needed."
            )

    # ---------------- PRODUCTION TAB ----------------

    def _createProductionTab(self):
        tab = Frame(self.tabs, setLayout=True)
        self.tabs.addTab(tab, "Production")

        row = 0
        Label(tab, "Include NMR restraints (DISANG):", grid=(row, 0))
        self.prodUseRst = RadioButtons(
            tab, texts=["No", "Yes"], selectedInd=0, grid=(row, 1), gridSpan=(1, 2)
        )
        row += 1

        Label(tab, "Local CPUs:", grid=(row, 0))
        self.prodNCpu = Entry(tab, text="8", grid=(row, 1))
        row += 1

        Label(tab, "Local GPUs (0/1):", grid=(row, 0))
        self.prodLocalGpus = Entry(tab, text="0", grid=(row, 1))
        row += 1

        Label(tab, "Steps (nstlim):", grid=(row, 0))
        self.prodSteps = Entry(tab, text="250000", grid=(row, 1))
        row += 1

        Label(tab, "Timestep (dt):", grid=(row, 0))
        self.prodDt = Entry(tab, text="0.002", grid=(row, 1))
        row += 1

        Label(tab, "Start restart (.rst7):", grid=(row, 0))
        self.prodStartRst = Entry(tab, text="", grid=(row, 1))
        Button(tab, "Browse", callback=self._browseProdStart, grid=(row, 2))
        row += 1

        Label(tab, "Output subdir:", grid=(row, 0))
        self.prodOut = Entry(tab, text="03_production", grid=(row, 1))
        row += 1

        Button(
            tab,
            "Run production locally",
            callback=self._runProduction,
            grid=(row, 0),
            gridSpan=(1, 3),
        )
        row += 1

        HLine(tab, grid=(row, 0), gridSpan=(1, 3))
        row += 1

        Label(tab, "HTCondor CPUs:", grid=(row, 0))
        self.prodCondorCpus = Entry(tab, text="16", grid=(row, 1))
        row += 1

        Label(tab, "HTCondor GPUs:", grid=(row, 0))
        self.prodCondorGpus = Entry(tab, text="0", grid=(row, 1))
        row += 1

        Label(tab, "HTCondor memory (GB):", grid=(row, 0))
        self.prodCondorMemGb = Entry(tab, text="32", grid=(row, 1))
        row += 1

        Label(tab, "HTCondor MaxRuntime (s):", grid=(row, 0))
        self.prodCondorMaxRt = Entry(tab, text="259200", grid=(row, 1))
        row += 1

        Button(
            tab,
            "Generate HTCondor scripts",
            callback=self._condorProdGenerate,
            grid=(row, 0),
        )
        Button(
            tab, "Submit to HTCondor", callback=self._condorProdSubmit, grid=(row, 1)
        )
        Button(tab, "Check job status", callback=self._condorProdStatus, grid=(row, 2))

    def _browseProdStart(self):
        f, _ = QFileDialog.getOpenFileName(
            self,
            "Select restart (rst7)",
            filter="Restart (*.rst7 *.rst *.nc);;All (*.*)",
        )
        if f:
            self.prodStartRst.setText(f)

    def _runProduction(self):
        try:
            wd = self._getWorkDir()
            if not wd:
                self._log("ERROR: Work dir is empty (set it in Input tab)")
                return

            # write helper script
            outsub = str(self.prodOut.get()).strip() or "03_production"
            script_path = os.path.join(wd, "production.sh")
            try:
                with open(script_path, "w") as fh:
                    fh.write(PROD_SH + "\n")
                os.chmod(script_path, 0o755)
            except Exception as e:
                self._log(f"ERROR: could not write production.sh: {e}")
                return

            # Local CPUs
            try:
                np = int(str(self.prodNCpu.get()).strip())
                if np < 1:
                    np = 1
            except Exception:
                np = 1

            # Local GPUs (0/1 only)
            try:
                ngpus = int(str(self.prodLocalGpus.get()).strip())
                if ngpus < 0:
                    ngpus = 0
                if ngpus > 1:
                    ngpus = 1
            except Exception:
                ngpus = 0

            # Steps / dt
            try:
                nstlim = int(str(self.prodSteps.get()).strip())
                if nstlim < 1:
                    nstlim = 250000
            except Exception:
                nstlim = 250000

            try:
                dt = float(str(self.prodDt.get()).strip())
                if dt <= 0.0:
                    dt = 0.002
            except Exception:
                dt = 0.002

            use_rst = "1" if self.prodUseRst.get() == 1 else "0"

            start_rst = str(self.prodStartRst.get()).strip()

            cmd = [
                "bash",
                script_path,
                "--workdir",
                wd,
                "--np",
                str(np),
                "--ngpus",
                str(ngpus),
                "--temp",
                "300.0",
                "--nstlim",
                str(nstlim),
                "--dt",
                str(dt),
                "--use_rst",
                use_rst,
                "--outdir",
                outsub,
            ]
            if start_rst:
                cmd += ["--start_rst", start_rst]

            self._runCommand(cmd, tag="PRODUCTION_LOCAL")

        except Exception as e:
            self._log(f"[PROD:CONDOR] ERROR in _runProduction: {e}")
            raise

    def _condor_prod_dir(self, wd):
        return os.path.join(wd, "03_production", "htcondor_production")

    def _condorProdGenerate(self):
        wd = self._getWorkDir()
        if not wd:
            self._log("ERROR: Work dir is empty (set it in Input tab)")
            return

        def _int(widget, default, minv=None):
            try:
                v = int(str(widget.get()).strip())
                if minv is not None and v < minv:
                    return default
                return v
            except Exception:
                return default

        req_cpus = _int(self.prodCondorCpus, 16, 1)
        req_gpus = _int(self.prodCondorGpus, 0, 0)
        mem_gb = _int(self.prodCondorMemGb, 32, 1)
        maxrt = _int(self.prodCondorMaxRt, 259200, 60)

        use_rst = "1" if self.prodUseRst.get() == 1 else "0"

        htdir = self._condor_prod_dir(wd)
        os.makedirs(htdir, exist_ok=True)
        os.makedirs(os.path.join(htdir, "OUT"), exist_ok=True)

        run_path = os.path.join(htdir, "run_pipeline.sh")
        prod_path = os.path.join(htdir, "production.sh")
        submit_path = os.path.join(htdir, "amber_pmemd_production.submit")

        # Stage preprocess inputs (relative paths preserved)
        prep_src = os.path.join(wd, "01_preprocess")
        prep_dst = os.path.join(htdir, "01_preprocess")
        os.makedirs(prep_dst, exist_ok=True)
        for fn in ("input_salt.prmtop", "input_salt.inpcrd", "input.RST"):
            src = os.path.join(prep_src, fn)
            if os.path.isfile(src):
                shutil.copy2(src, os.path.join(prep_dst, fn))

        # Stage starting restart for production
        start_rst = str(self.prodStartRst.get()).strip()
        if not start_rst:
            cand1 = os.path.join(wd, "02_minrelax", "OUT", "9md.rst7")
            cand2 = os.path.join(wd, "02_minrelax", "9md.rst7")
            start_rst = cand1 if os.path.isfile(cand1) else cand2
        if not os.path.isfile(start_rst):
            self._log(f"ERROR: starting restart not found for production: {start_rst}")
            return
        shutil.copy2(start_rst, os.path.join(htdir, "start.rst7"))

        nstlim = str(self.prodSteps.get()).strip() or "250000"
        dt = str(self.prodDt.get()).strip() or "0.002"

        wrapper = f"""#!/bin/bash
set -e

# HTCondor passes NPROCS/NGPUS via 'environment' in the submit file.
NPROCS="${{NPROCS:-{req_cpus}}}"
NGPUS="${{NGPUS:-{req_gpus}}}"

# MPI engines require >=2 ranks
if [[ "$NGPUS" -ge 2 ]] && [[ "$NPROCS" -lt 2 ]]; then NPROCS=2; fi
if [[ "$NGPUS" -eq 0 ]] && [[ "$NPROCS" -lt 2 ]]; then NPROCS=2; fi

OUT="OUT"
mkdir -p "$OUT"
echo "Job started: $(date) on $(hostname)" > "$OUT/_job_meta.txt"

if [ -r /usr/software/amber/amber.sh ]; then
  source /usr/software/amber/amber.sh || true
fi

bash ./production.sh \
  --workdir "." \
  --np "$NPROCS" \
  --ngpus "$NGPUS" \
  --temp "300.0" \
  --nstlim "{nstlim}" \
  --dt "{dt}" \
  --use_rst "{use_rst}" \
  --rstfile "./01_preprocess/input.RST" \
  --outdir "OUT" \
  --start_rst "./start.rst7"

echo "Job finished: $(date)" >> "$OUT/_job_meta.txt"
"""
        with open(run_path, "w") as fh:
            fh.write(wrapper)
        os.chmod(run_path, 0o755)

        with open(prod_path, "w") as fh:
            fh.write(PROD_SH + "\n")
        os.chmod(prod_path, 0o755)

        submit_lines = [
            "universe                = vanilla",
            "executable              = /bin/bash",
            "arguments               = run_pipeline.sh",
            "",
            f"request_cpus            = {req_cpus}",
            f"request_gpus            = {req_gpus}",
            f"request_memory          = {mem_gb} GB",
            f"+MaxRuntime             = {maxrt}",
            "",
            'environment             = "NPROCS=$(RequestCpus) NGPUS=$(RequestGpus)"',
            "getenv                  = True",
            "",
            "transfer_executable     = false",
            "should_transfer_files   = YES",
            "when_to_transfer_output = ON_EXIT",
            "",
            "transfer_input_files    = run_pipeline.sh, production.sh, start.rst7, 01_preprocess/input_salt.prmtop, 01_preprocess/input_salt.inpcrd, 01_preprocess/input.RST",
            "transfer_output_files   = OUT",
            "",
            'requirements            = (AMBER =!= undefined) && (AMBERTOOLS =!= undefined) && regexp("nmrbox.org", Machine)',
            "log                     = condor.log",
            "output                  = condor.out",
            "error                   = condor.err",
            "",
            "queue 1",
        ]
        with open(submit_path, "w") as fh:
            fh.write("\n".join(submit_lines) + "\n")

        self._log("[CONDOR:PRODUCTION] Generated HTCondor scripts")
        self._log(f"[CONDOR:PRODUCTION] dir={htdir}")

    def _condorProdSubmit(self):
        wd = self._getWorkDir()
        if not wd:
            self._log("ERROR: Work dir is empty (set it in Input tab)")
            return
        htdir = self._condor_prod_dir(wd)
        submit_path = os.path.join(htdir, "amber_pmemd_production.submit")
        if not os.path.isfile(submit_path):
            self._log("ERROR: submit file not found. Generate HTCondor scripts first.")
            return
        self._runCommand(
            ["condor_submit", submit_path], tag="CONDOR_PRODUCTION", cwd=htdir
        )

    def _condorProdStatus(self):
        self._run("CONDOR_Q", "condor_q -submitter $USER")

    def _createAnalysisTab(self):
        tab = Frame(setLayout=True)
        self.tabs.addTab(tab, "Analysis")

        row = 0
        Label(tab, "Topology (prmtop):", grid=(row, 0))
        self.anPrmtop = Entry(tab, grid=(row, 1))
        Button(tab, "Browse", grid=(row, 2), callback=self._browseAnPrmtop)
        row += 1

        Label(tab, "Trajectory (nc):", grid=(row, 0))
        self.anTraj = Entry(tab, grid=(row, 1))
        Button(tab, "Browse", grid=(row, 2), callback=self._browseAnTraj)
        row += 1

        Label(tab, "MD output (prod.out):", grid=(row, 0))
        self.anMdout = Entry(tab, grid=(row, 1))
        Button(tab, "Browse", grid=(row, 2), callback=self._browseAnMdout)
        row += 1

        Label(tab, "Atom selection:", grid=(row, 0))
        self.anSel = PulldownList(
            tab,
            texts=[
                "Backbone (C,CA,N,O)",
                "C-alpha (CA)",
                "Protein heavy atoms",
                "All protein atoms",
            ],
            grid=(row, 1),
        )
        row += 1

        Label(tab, "Output folder name:", grid=(row, 0))
        self.anOut = Entry(tab, text="04_analysis", grid=(row, 1))
        row += 1

        self.doRmsd = CheckBox(tab, text="RMSD", checked=True, grid=(row, 0))
        self.doRmsf = CheckBox(tab, text="RMSF", checked=True, grid=(row, 1))
        row += 1
        self.doViol = CheckBox(
            tab, text="Restraint violations (from mdout)", checked=True, grid=(row, 0)
        )
        self.doPlots = CheckBox(
            tab, text="Export plots (PNG)", checked=True, grid=(row, 1)
        )
        row += 1

        Button(tab, "Run selected analyses", grid=(row, 1), callback=self._runAnalysis)
        row += 1
        Button(
            tab,
            "View trajectory in PyMOL",
            grid=(row, 1),
            callback=self._viewTrajectoryInPymol,
        )

        self._setAnalysisDefaults()

    def _setAnalysisDefaults(self):
        wd = self.workDir.get()
        if not wd:
            return
        prmtop = os.path.join(wd, "01_preprocess", "input_salt.prmtop")
        traj = os.path.join(wd, "03_production", "prod.nc")
        mdout = os.path.join(wd, "03_production", "prod.out")
        if os.path.exists(prmtop) and not self.anPrmtop.get():
            self.anPrmtop.setText(prmtop)
        if os.path.exists(traj) and not self.anTraj.get():
            self.anTraj.setText(traj)
        if os.path.exists(mdout) and not self.anMdout.get():
            self.anMdout.setText(mdout)

    def _browseAnPrmtop(self):
        f, _ = QFileDialog.getOpenFileName(
            self, "Select topology", filter="Topology (*.prmtop *.parm7);;All (*.*)"
        )
        if f:
            self.anPrmtop.setText(f)

    def _browseAnTraj(self):
        f, _ = QFileDialog.getOpenFileName(
            self,
            "Select trajectory",
            filter="Trajectory (*.nc *.mdcrd *.crd *.dcd *.xtc);;All (*.*)",
        )
        if f:
            self.anTraj.setText(f)

    def _browseAnMdout(self):
        f, _ = QFileDialog.getOpenFileName(
            self, "Select mdout", filter="Output (*.out *.mdout *.log);;All (*.*)"
        )
        if f:
            self.anMdout.setText(f)

    def _cpptraj_mask_from_choice(self, ind):
        if ind == 0:
            return ":1-999@C,CA,N,O"
        elif ind == 1:
            return ":1-999@CA"
        elif ind == 2:
            return ":1-999 & !@H="
        else:
            return ":1-999"

    def _runAnalysis(self):
        # Robust analysis runner with defensive logging
        self._log("[ANALYSIS] Button clicked")
        try:
            # --- Resolve workdir (allow analysis-only use) ---
            wd = (self.workDir.get() or "").strip()
            if not wd:
                self._log(
                    "[ANALYSIS] Workdir is empty; trying to infer from provided file paths"
                )
                cand_paths = []
                for w in (
                    getattr(self, "anPrmtop", None),
                    getattr(self, "anTraj", None),
                    getattr(self, "anMdout", None),
                ):
                    try:
                        p = (w.get() or "").strip() if w is not None else ""
                    except Exception:
                        p = ""
                    if p:
                        cand_paths.append(os.path.abspath(os.path.expanduser(p)))
                if cand_paths:
                    try:
                        common = os.path.commonpath(cand_paths)
                    except Exception:
                        common = os.path.dirname(cand_paths[0])
                    if os.path.isfile(common):
                        common = os.path.dirname(common)
                    wd = common
                    try:
                        self.workDir.setText(wd)
                    except Exception:
                        pass
                    self._log(f"[ANALYSIS] Inferred workdir: {wd}")
                else:
                    showWarning(
                        "Missing input",
                        "Set Working directory or provide prmtop/trajectory paths so I can infer it.",
                    )
                    self._log(
                        "[ANALYSIS] No workdir and no file paths to infer from; aborting"
                    )
                    return

            # Create output folder early so user sees something even if later checks fail
            outdir = (
                (self.anOut.get() or "04_analysis").strip()
                if hasattr(self, "anOut")
                else "04_analysis"
            )
            abs_outdir = os.path.join(wd, outdir)
            os.makedirs(abs_outdir, exist_ok=True)
            self._log(f"[ANALYSIS] Output directory: {abs_outdir}")

            # Refresh auto-defaults (workdir may have been set after tab creation)
            if hasattr(self, "_setAnalysisDefaults"):
                self._setAnalysisDefaults()

            # --- Collect inputs ---
            def _guess(p):
                return p if p and os.path.exists(p) else ""

            prmtop = (
                (self.anPrmtop.get() or "").strip() if hasattr(self, "anPrmtop") else ""
            )
            traj = (self.anTraj.get() or "").strip() if hasattr(self, "anTraj") else ""
            mdout = (
                (self.anMdout.get() or "").strip() if hasattr(self, "anMdout") else ""
            )

            # Common defaults if not specified
            if not prmtop:
                prmtop = _guess(os.path.join(wd, "01_preprocess", "input_salt.prmtop"))
            if not traj:
                traj = _guess(os.path.join(wd, "03_production", "prod.nc"))
            if not mdout:
                mdout = _guess(os.path.join(wd, "03_production", "prod.out"))

            # Write back guessed paths so user can see them
            if (
                hasattr(self, "anPrmtop")
                and prmtop
                and not (self.anPrmtop.get() or "").strip()
            ):
                self.anPrmtop.setText(prmtop)
            if (
                hasattr(self, "anTraj")
                and traj
                and not (self.anTraj.get() or "").strip()
            ):
                self.anTraj.setText(traj)
            if (
                hasattr(self, "anMdout")
                and mdout
                and not (self.anMdout.get() or "").strip()
            ):
                self.anMdout.setText(mdout)

            self._log(f"[ANALYSIS] prmtop: {prmtop}")
            self._log(f"[ANALYSIS] traj  : {traj}")
            self._log(f"[ANALYSIS] mdout : {mdout if mdout else '(not set)'}")

            # Validate files
            if not prmtop or not os.path.exists(prmtop):
                showWarning("Missing file", "Topology (prmtop) not set or not found")
                self._log(f"[ANALYSIS] ERROR: prmtop missing/not found: {prmtop}")
                return
            if not traj or not os.path.exists(traj):
                showWarning("Missing file", "Trajectory (.nc) not set or not found")
                self._log(f"[ANALYSIS] ERROR: traj missing/not found: {traj}")
                return

            # --- Determine selection mask ---
            sel_mask = ":1-999@C,CA,N,O"
            if hasattr(self, "anSel"):
                try:
                    sel_mask = self._cpptraj_mask_from_choice(self._idx(self.anSel))
                except Exception:
                    pass
            self._log(f"[ANALYSIS] cpptraj selection mask: {sel_mask}")

            # --- Ensure scripts exist ---
            analysis_sh = self._ensure_script(
                wd, "analysis_cpptraj.sh", ANALYSIS_CPPTRAJ_SH
            )
            parser_py = self._ensure_script(
                wd, "parse_violations.py", PARSER_VIOLATIONS_PY
            )
            plot_py = self._ensure_script(wd, "plot_analysis.py", PLOT_ANALYSIS_PY)

            # --- Determine what to run (defensive against missing widgets) ---
            def _checked(name):
                w = getattr(self, name, None)
                try:
                    return bool(w and w.isChecked())
                except Exception:
                    return False

            do_rmsd = _checked("doRmsd")
            do_rmsf = _checked("doRmsf")
            do_viol = _checked("doViol")
            do_plots = _checked("doPlots")

            self._log(
                f"[ANALYSIS] selected -> RMSD:{do_rmsd} RMSF:{do_rmsf} Violations:{do_viol} Plots:{do_plots}"
            )

            # --- Build commands ---
            cmds = []

            # cpptraj needs AMBER environment
            if do_rmsd or do_rmsf:
                cmds.append(
                    f"source '{AMBER_SH}' && '{analysis_sh}' --prmtop '{prmtop}' --traj '{traj}' --sel '{sel_mask}' --outdir '{outdir}'"
                )

            # Python-only steps; activate venv if available
            venv = ""
            if hasattr(self, "venvPath"):
                try:
                    venv = (self.venvPath.get() or "").strip()
                except Exception:
                    venv = ""
            if not venv:
                venv = DEFAULT_VENV

            venv_activate = ""
            if venv and os.path.exists(
                os.path.join(os.path.expanduser(venv), "bin", "activate")
            ):
                venv_activate = f"source '{os.path.expanduser(venv)}/bin/activate' && "

            if do_viol:
                if mdout and os.path.exists(mdout):
                    cmds.append(
                        f"{venv_activate}python3 '{parser_py}' --mdout '{mdout}' --out '{os.path.join(abs_outdir, 'violations.csv')}'"
                    )
                else:
                    self._log(
                        "[ANALYSIS] NOTE: mdout not found; skipping violation parsing"
                    )

            if do_plots:
                # Use system python for plotting: matplotlib may be available in the system environment
                # but not installed in the pdb2pqr venv.
                cmds.append(f"python3 '{plot_py}' --analysis_dir '{abs_outdir}'")

            if not cmds:
                self._log(
                    "[ANALYSIS] Nothing selected to run. Tick RMSD/RMSF/Violations/Plots."
                )
                showWarning(
                    "Nothing selected",
                    "Tick at least one analysis option (RMSD/RMSF/Violations/Plots).",
                )
                return

            cmd = "set -e; " + " ; ".join(cmds)
            self._log("[RUN:ANALYSIS] " + cmd)
            self._run("ANALYSIS", cmd, wd)

        except Exception as e:
            import traceback

            tb = traceback.format_exc()
            self._log("[ANALYSIS] EXCEPTION: " + str(e))
            self._log(tb)
            showWarning(
                "Analysis error",
                f"Analysis crashed with:\n{e}\n\nSee Command log for traceback.",
            )

    def _viewTrajectoryInPymol(self):
        """Launch PyMOL to visualise the provided AMBER trajectory.

        Stable NMRbox/VNC recipe (confirmed by user):
          - unset PYTHONHOME/PYTHONPATH/PYTHONUSERBASE and set PYTHONNOUSERSITE=1
          - force software OpenGL (LIBGL_ALWAYS_SOFTWARE=1) to avoid VirtualGL PBO errors
          - start PyMOL with -k to skip user rc/plugins (avoids bus errors)

        Additionally, auto-generate a matching reference PDB (frame0.pdb) using cpptraj
        so the NetCDF trajectory loads reliably.
        """
        self._log("[PYMOL] Button clicked")
        try:
            wd = (self.workDir.get() or "").strip()
            prmtop = (
                (self.anPrmtop.get() or "").strip() if hasattr(self, "anPrmtop") else ""
            )
            traj = (self.anTraj.get() or "").strip() if hasattr(self, "anTraj") else ""

            # Infer workdir if not set
            if not wd:
                for p in (prmtop, traj):
                    if p:
                        wd = pathlib.Path(p).resolve().parent.as_posix()
                        if pathlib.Path(wd).name in (
                            "01_preprocess",
                            "02_minrelax",
                            "03_production",
                            "04_analysis",
                            "05_analysis",
                            "scripts",
                        ):
                            wd = pathlib.Path(wd).parent.as_posix()
                        break

            if not wd:
                showWarning(
                    "PyMOL",
                    "Working directory is not set and cannot be inferred.\nPlease set Workdir or provide valid prmtop/traj paths.",
                )
                return

            if not prmtop or not os.path.isfile(prmtop):
                showWarning("PyMOL", f"Topology file not found:\n{prmtop}")
                return
            if not traj or not os.path.isfile(traj):
                showWarning("PyMOL", f"Trajectory file not found:\n{traj}")
                return

            prmtop_abs = os.path.abspath(prmtop)
            traj_abs = os.path.abspath(traj)

            # Analysis dir for frame0.pdb
            outdir = os.path.join(wd, "05_analysis")
            os.makedirs(outdir, exist_ok=True)
            frame0_pdb = os.path.join(outdir, "frame0.pdb")

            scripts_dir = os.path.join(wd, "scripts")
            os.makedirs(scripts_dir, exist_ok=True)

            # 1) Generate frame0.pdb with cpptraj
            cpp_in = os.path.join(scripts_dir, "cpptraj_frame0.in")
            with open(cpp_in, "w") as fh:
                fh.write(f"trajin {traj_abs} 1 1\n")
                fh.write(f"trajout {frame0_pdb} pdb\n")
                fh.write("go\n")

            cpp_cmd = "set -e; source '/usr/software/amber/amber.sh' && cpptraj -p '{p}' -i '{i}'".format(
                p=prmtop_abs, i=cpp_in
            )
            self._log(f"[PYMOL] Generating frame0 PDB with cpptraj -> {frame0_pdb}")
            self._log("[RUN:CPPTRAJ] " + cpp_cmd)

            # Run cpptraj synchronously (QProcess is async and can trigger a premature existence check)
            try:
                proc = subprocess.run(
                    [SHELL, "-lc", cpp_cmd], cwd=wd, text=True, capture_output=True
                )
                if proc.stdout:
                    for line in proc.stdout.splitlines():
                        self._log(line)
                if proc.stderr:
                    for line in proc.stderr.splitlines():
                        self._log(line)
                if proc.returncode != 0:
                    showWarning(
                        "PyMOL",
                        f"cpptraj failed (exit code {proc.returncode}).\nSee Command log for details.",
                    )
                    return
            except Exception as e:
                showWarning("PyMOL", f"cpptraj execution failed:\n{e}")
                return

            if not os.path.isfile(frame0_pdb) or os.path.getsize(frame0_pdb) < 100:
                showWarning(
                    "PyMOL",
                    f"frame0.pdb was not generated correctly:\n{frame0_pdb}\n\nCheck Command log for cpptraj errors.",
                )
                return

            # 2) Write PyMOL script that loads frame0.pdb then the NetCDF trajectory
            pml_path = os.path.join(scripts_dir, "view_traj.pml")
            pml = f"""reinitialize
set auto_zoom, 0
load {frame0_pdb}, mol
load_traj {traj_abs}, mol, state=1
hide everything, mol
show cartoon, mol
set cartoon_sampling, 5
bg_color white
zoom mol
"""
            with open(pml_path, "w") as fh:
                fh.write(pml)

            # 3) Launch PyMOL detached; capture log
            log_path = os.path.join(scripts_dir, "pymol_launch.log")
            launch = (
                "unset PYTHONHOME PYTHONPATH PYTHONUSERBASE; "
                "export PYTHONNOUSERSITE=1; "
                "export LIBGL_ALWAYS_SOFTWARE=1; "
                f"nohup /usr/software/pymol/pymol -qk '{pml_path}' > '{log_path}' 2>&1 & "
                "echo '[PYMOL] PID='$!; sleep 3; "
                f"echo '[PYMOL] last log lines:'; tail -n 120 '{log_path}' || true"
            )

            self._log(f"[PYMOL] frame0.pdb: {frame0_pdb}")
            self._log(f"[PYMOL] PML written: {pml_path}")
            self._log(f"[PYMOL] Log file: {log_path}")
            self._log("[RUN:PYMOL] " + launch)
            self._run("PYMOL", launch, wd)

        except Exception as e:
            import traceback

            tb = traceback.format_exc()
            self._log("[PYMOL] EXCEPTION: " + str(e))
            self._log(tb)
            showWarning(
                "PyMOL error",
                f"Could not start PyMOL:\n{e}\n\nSee Command log for details.",
            )


try:
    mw = mainWindow
except NameError:
    from ccpn.framework.Application import getApplication

    mw = getApplication().mainWindow

popup = MdWorkflowDialog(parent=mw, mainWindow=mw)
popup.show()
popup.raise_()

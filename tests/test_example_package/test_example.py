# tests/test_NMR_MD/test_example.py
import pytest


class TestNMRMD:
    """
    Class-based placeholder tests for NMR_MD.

    - Demonstrates using a test class with pytest.
    - Can easily be extended with more methods for real tests.
    """

    def test_mdworkflowdialog_imports(self):
        """
        Basic smoke test to ensure MdWorkflowDialog can be imported.
        """
        pytest.importorskip("ccpn")  # skip cleanly if dependency missing
        from NMR_MD.macro95_final import MdWorkflowDialog

        assert MdWorkflowDialog is not None

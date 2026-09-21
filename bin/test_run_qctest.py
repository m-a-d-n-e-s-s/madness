#!/usr/bin/env python3
"""Unit tests for run_qctest.compare: the tol / rtol / max / allow_zero semantics
every qctest check.json relies on.

Run by hand:   python3 bin/test_run_qctest.py
Via ctest:     madness/test/qc/run_qctest_selftest/run   (labels qctest;short)
"""
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_qctest  # noqa: E402


class CompareTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        d = Path(self.tmp.name)
        self.out, self.ref = d / "out.json", d / "ref.json"

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, out, ref):
        self.out.write_text(json.dumps(out))
        self.ref.write_text(json.dumps(ref))

    def compare(self, checks):
        return run_qctest.compare(self.out, self.ref, checks)

    def test_tol_passes_within_and_fails_outside(self):
        self.write({"e": 1.0005}, {"e": 1.0})
        self.assertTrue(self.compare([{"key": ["e"], "tol": 1e-3}]))
        self.assertFalse(self.compare([{"key": ["e"], "tol": 1e-4}]))

    def test_rtol_is_relative_to_the_reference(self):
        self.write({"b": 101.0}, {"b": 100.0})
        self.assertTrue(self.compare([{"key": ["b"], "rtol": 0.02}]))
        self.assertFalse(self.compare([{"key": ["b"], "rtol": 0.005}]))

    def test_rtol_rejects_a_zero_reference_unless_allow_zero(self):
        self.write({"b": 0.0}, {"b": 0.0})
        self.assertFalse(self.compare([{"key": ["b"], "rtol": 0.01}]))
        self.assertTrue(self.compare([{"key": ["b"], "rtol": 0.01, "allow_zero": True}]))

    def test_max_bounds_the_output_only(self):
        self.write({"it": 9}, {"it": 7})  # an exact int tol would fail; max must not
        self.assertTrue(self.compare([{"key": ["it"], "max": 12}]))
        self.assertFalse(self.compare([{"key": ["it"], "max": 8}]))

    def test_max_rejects_a_non_numeric_value(self):
        self.write({"s": "converged"}, {"s": "converged"})
        self.assertFalse(self.compare([{"key": ["s"], "max": 1}]))

    def test_max_and_tol_both_apply(self):
        self.write({"it": 9}, {"it": 9})
        self.assertTrue(self.compare([{"key": ["it"], "max": 12, "tol": 0}]))
        self.write({"it": 9}, {"it": 8})
        self.assertFalse(self.compare([{"key": ["it"], "max": 12, "tol": 0}]))

    def test_missing_key_fails_for_every_kind(self):
        self.write({"a": 1.0}, {"a": 1.0})
        self.assertFalse(self.compare([{"key": ["zz"], "max": 1}]))
        self.assertFalse(self.compare([{"key": ["zz"], "rtol": 0.1}]))
        self.assertFalse(self.compare([{"key": ["zz"], "tol": 0.1}]))

    def test_strings_and_bools_still_compare_exactly(self):
        self.write({"s": "converged", "c": True}, {"s": "converged", "c": True})
        self.assertTrue(self.compare([{"key": ["s"], "tol": 0}, {"key": ["c"], "tol": 0}]))
        self.write({"s": "unconverged", "c": False}, {"s": "converged", "c": True})
        self.assertFalse(self.compare([{"key": ["s"], "tol": 0}]))
        self.assertFalse(self.compare([{"key": ["c"], "tol": 0}]))

    def test_zero_reference_with_tol_needs_allow_zero(self):
        self.write({"w": 0.0}, {"w": 0.0})
        self.assertFalse(self.compare([{"key": ["w"], "tol": 0}]))
        self.assertTrue(self.compare([{"key": ["w"], "tol": 0, "allow_zero": True}]))


if __name__ == "__main__":
    unittest.main()

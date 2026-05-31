"""Tests for shared AFIR helper behavior."""

import unittest

from pyar.biases.afir import resolve_gamma


class AfirUtilsTests(unittest.TestCase):
    def test_resolve_gamma_uses_fallback_for_missing_value(self):
        self.assertEqual(resolve_gamma(None), 100.0)

    def test_resolve_gamma_parses_numeric_strings(self):
        self.assertEqual(resolve_gamma("37.5"), 37.5)

    def test_resolve_gamma_rejects_invalid_value(self):
        with self.assertRaisesRegex(ValueError, "Invalid AFIR gamma"):
            resolve_gamma("not-a-number", fallback=25.0)

    def test_resolve_gamma_rejects_nonfinite_value(self):
        with self.assertRaisesRegex(ValueError, "Invalid AFIR gamma"):
            resolve_gamma("nan")


if __name__ == "__main__":
    unittest.main()

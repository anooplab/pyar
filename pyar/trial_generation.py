"""Compatibility shim for the legacy ``pyar.trial_generation`` module."""

if __name__ == "__main__":
    import runpy

    runpy.run_module("pyar.sampling.trial_generator", run_name="__main__")
else:
    from pyar.sampling.trial_generator import *  # noqa: F401,F403

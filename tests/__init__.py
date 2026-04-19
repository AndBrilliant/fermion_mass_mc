"""tests/ — v4 MC test modules.

Running order for a production run:
    python -m tests.validate_reproducibility       # smoke test, few seconds
    python -m tests.test_a_random_spectrum         # Null A, ~minutes
    python -m tests.test_b_bootstrap               # Null B, ~minute
    python -m tests.test_c_algorithm_space         # Algorithm check, ~minute
    python -m tests.test_d_temporal_convergence    # Temporal, ~seconds

Or chain them:
    python -m tests.validate_reproducibility && \\
    python -m tests.test_a_random_spectrum && \\
    python -m tests.test_b_bootstrap && \\
    python -m tests.test_c_algorithm_space && \\
    python -m tests.test_d_temporal_convergence
"""

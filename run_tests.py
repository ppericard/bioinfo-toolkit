#!/usr/bin/env python3
"""
Script to run tests for bioinfo-toolkit.

This script provides a convenient way to run different types of tests
for the bioinfo-toolkit repository.
"""

import os
import sys
import subprocess
import argparse


def run_unit_tests():
    """Run unit tests."""
    print("Running unit tests...")
    subprocess.run(["pytest", "tests/unit", "-v"], check=True)


def run_integration_tests():
    """Run integration tests."""
    print("Running integration tests...")
    subprocess.run(["pytest", "tests/integration", "-v"], check=True)


def run_functional_tests():
    """Run functional tests."""
    print("Running functional tests...")
    subprocess.run(["pytest", "tests/functional", "-v"], check=True)


def run_performance_tests():
    """Run performance tests."""
    print("Running performance tests...")
    subprocess.run(["pytest", "tests/performance", "-v"], check=True)


def run_all_tests(with_coverage=False):
    """Run all tests."""
    print("Running all tests...")
    if with_coverage:
        subprocess.run(
            ["pytest", "--cov=src", "--cov-report=term", "--cov-report=html"],
            check=True,
        )
    else:
        subprocess.run(["pytest"], check=True)


def run_linting():
    """Run code linting."""
    print("Running code linting...")
    subprocess.run(
        ["flake8", "--count", "--max-complexity=10", "--max-line-length=100", "src", "tests"],
        check=True,
    )


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Run tests for bioinfo-toolkit")
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--unit", action="store_true", help="Run unit tests"
    )
    group.add_argument(
        "--integration", action="store_true", help="Run integration tests"
    )
    group.add_argument(
        "--functional", action="store_true", help="Run functional tests"
    )
    group.add_argument(
        "--performance", action="store_true", help="Run performance tests"
    )
    group.add_argument(
        "--lint", action="store_true", help="Run code linting"
    )
    parser.add_argument(
        "--coverage", action="store_true", help="Generate coverage report"
    )
    
    args = parser.parse_args()
    
    if args.unit:
        run_unit_tests()
    elif args.integration:
        run_integration_tests()
    elif args.functional:
        run_functional_tests()
    elif args.performance:
        run_performance_tests()
    elif args.lint:
        run_linting()
    else:
        run_all_tests(with_coverage=args.coverage)
    
    return 0


if __name__ == "__main__":
    sys.exit(main()) 
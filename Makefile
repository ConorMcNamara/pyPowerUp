.PHONY: help install install-dev test test-cov lint format type-check clean build publish pre-commit requirements

help:
	@echo "Available commands:"
	@echo "  make install        Install package in production mode"
	@echo "  make install-dev    Install package in development mode with dev dependencies"
	@echo "  make test           Run tests"
	@echo "  make test-cov       Run tests with coverage report"
	@echo "  make lint           Run linter (ruff)"
	@echo "  make format         Format code with ruff"
	@echo "  make type-check     Run type checker (mypy)"
	@echo "  make clean          Remove build artifacts and cache files"
	@echo "  make build          Build distribution packages"
	@echo "  make publish        Publish to PyPI (requires credentials)"
	@echo "  make pre-commit     Install and run pre-commit hooks"
	@echo "  make requirements   Generate requirements.txt from pyproject.toml (legacy support)"
	@echo "  make all            Run format, lint, type-check, and test"

install:
	pip install -e .

install-dev:
	pip install -e ".[dev]"
	pre-commit install

test:
	pytest -v

test-cov:
	pytest --cov=pyPowerUp --cov-report=html --cov-report=term-missing

lint:
	ruff check .

format:
	ruff format .
	ruff check --fix .

type-check:
	mypy pyPowerUp --ignore-missing-imports

clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info
	rm -rf .pytest_cache/
	rm -rf .mypy_cache/
	rm -rf .ruff_cache/
	rm -rf htmlcov/
	rm -rf .coverage
	rm -rf coverage.xml
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name '*.pyc' -delete
	find . -type f -name '*.pyo' -delete
	find . -type f -name '*~' -delete

build: clean
	python -m build

publish: build
	python -m twine upload dist/*

pre-commit:
	pre-commit install
	pre-commit run --all-files

requirements:
	@echo "# This file is auto-generated from pyproject.toml" > requirements.txt
	@echo "# To regenerate: make requirements" >> requirements.txt
	@echo "# For modern Python projects, install directly from pyproject.toml:" >> requirements.txt
	@echo "#   pip install -e ." >> requirements.txt
	@echo "" >> requirements.txt
	@python3 -c "import tomllib; f = open('pyproject.toml', 'rb'); data = tomllib.load(f); deps = data['project']['dependencies']; print('\n'.join(deps))" >> requirements.txt
	@echo "Generated requirements.txt from pyproject.toml"

all: format lint type-check test
	@echo "All checks passed!"

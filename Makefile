.PHONY: help install install-dev test test-cov lint format type-check clean build publish pre-commit requirements

help:
	@echo "Available commands:"
	@echo "  make install        Install package in production mode"
	@echo "  make install-dev    Install package in development mode with dev dependencies"
	@echo "  make test           Run tests"
	@echo "  make test-cov       Run tests with coverage report"
	@echo "  make lint           Run linter (ruff)"
	@echo "  make format         Format code with ruff"
	@echo "  make type-check     Run type checker (zuban)"
	@echo "  make clean          Remove build artifacts and cache files"
	@echo "  make build          Build distribution packages"
	@echo "  make publish        Publish to PyPI (requires credentials)"
	@echo "  make pre-commit     Install and run pre-commit hooks"
	@echo "  make requirements   Export requirements.txt from uv.lock (legacy support)"
	@echo "  make all            Run format, lint, type-check, and test"

install:
	uv sync

install-dev:
	uv sync --extra dev
	uv run pre-commit install

test:
	uv run pytest -v

test-cov:
	uv run pytest --cov=pyPowerUp --cov-report=html --cov-report=term-missing

lint:
	uv run ruff check .

format:
	uv run ruff format .
	uv run ruff check --fix .

type-check:
	uv run zuban check pyPowerUp

clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info
	rm -rf .pytest_cache/
	rm -rf .zuban_cache/
	rm -rf .ruff_cache/
	rm -rf htmlcov/
	rm -rf .coverage
	rm -rf coverage.xml
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name '*.pyc' -delete
	find . -type f -name '*.pyo' -delete
	find . -type f -name '*~' -delete

build: clean
	uv build

publish: build
	uv publish

pre-commit:
	uv run pre-commit install
	uv run pre-commit run --all-files

requirements:
	uv export --no-dev --no-emit-project --format requirements-txt -o requirements.txt
	@echo "Exported requirements.txt from uv.lock"

all: format lint type-check test
	@echo "All checks passed!"

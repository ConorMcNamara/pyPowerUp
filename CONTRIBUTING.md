# Contributing to pyPowerUp

Thank you for your interest in contributing to pyPowerUp! This document provides guidelines and instructions for contributing to the project.

## Development Setup

### Prerequisites

- Python 3.13 or 3.14
- Git

### Setting Up Your Development Environment

1. Fork the repository on GitHub
2. Clone your fork locally:
   ```bash
   git clone https://github.com/YOUR-USERNAME/pyPowerUp.git
   cd pyPowerUp
   ```

3. Install development dependencies:
   ```bash
   pip install -e ".[dev]"
   ```

4. Install pre-commit hooks:
   ```bash
   pre-commit install
   ```

## Development Workflow

### Code Style

This project uses modern Python tooling for code quality:

- **Ruff** for linting and formatting
- **mypy** for type checking
- **pytest** for testing

All code must:
- Follow PEP 8 style guidelines (enforced by Ruff)
- Include type hints for all functions and methods
- Pass all linting and type checking
- Include NumPy-style docstrings

### Making Changes

1. Create a new branch for your feature or bugfix:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. Make your changes, ensuring:
   - Code is properly formatted: `make format`
   - Code passes linting: `make lint`
   - Type hints are correct: `make type-check`
   - Tests pass: `make test`

3. Add tests for new functionality in the `test/` directory

4. Update documentation as needed

### Managing Dependencies

This project uses `pyproject.toml` exclusively for dependency management (no `requirements.txt`). To add a new dependency:

1. Add it to the appropriate section in `pyproject.toml`:
   - Production dependencies → `dependencies`
   - Development dependencies → `project.optional-dependencies.dev`
   - Test dependencies → `project.optional-dependencies.test`

2. Reinstall the package:
   ```bash
   pip install -e ".[dev]"
   ```

Example:
```toml
[project]
dependencies = [
    "scipy>=1.15.2,<2.0.0",
    "numpy>=1.26.0,<3.0.0",
    "new-package>=1.0.0",  # Add here
]
```

### Running Tests

```bash
# Run all tests
make test

# Run tests with coverage
make test-cov

# Run all quality checks
make all
```

### Committing Your Changes

1. Pre-commit hooks will automatically run when you commit
2. Write clear, descriptive commit messages
3. Reference any relevant issues in your commit messages

Example:
```bash
git commit -m "Add power calculation for new design type

Implements power_xyz123 function for XYZ experimental design.
Closes #42"
```

## Pull Request Process

1. Push your branch to your fork:
   ```bash
   git push origin feature/your-feature-name
   ```

2. Create a Pull Request on GitHub with:
   - Clear description of the changes
   - Reference to any related issues
   - Examples of new functionality (if applicable)

3. Ensure all CI checks pass

4. Wait for review from maintainers

5. Address any requested changes

## Code Review

All submissions require review. We use GitHub pull requests for this purpose. Reviewers will look for:

- Code quality and style compliance
- Test coverage
- Documentation completeness
- Compatibility with existing code
- Performance considerations

## Testing Guidelines

- Write tests for all new functions
- Maintain or improve code coverage
- Test edge cases and error conditions
- Use descriptive test names

Example test structure:
```python
def test_power_calculation_basic():
    """Test basic power calculation with typical inputs."""
    result = power_bcra3f2(
        effect_size=0.145,
        rho2=0.10,
        n=20,
        J=44,
        K=5,
        alpha=0.05
    )
    assert 0.80 <= result <= 0.81
```

## Documentation

- Update docstrings for any changed functions
- Follow NumPy docstring style
- Include examples in docstrings when helpful
- Update README.md if adding major features

## Reporting Bugs

When reporting bugs, please include:

1. Python version
2. pyPowerUp version
3. Operating system
4. Minimal code example to reproduce the issue
5. Expected vs. actual behavior
6. Any error messages or stack traces

## Requesting Features

When requesting features:

1. Check if a similar feature request exists
2. Clearly describe the feature and its use case
3. Explain why it would be beneficial to the project
4. Consider offering to implement it yourself

## Code of Conduct

Be respectful and professional in all interactions. We aim to maintain a welcoming and inclusive environment for all contributors.

## Questions?

If you have questions, feel free to:
- Open an issue for discussion
- Contact the maintainers

## License

By contributing to pyPowerUp, you agree that your contributions will be licensed under the MIT License.

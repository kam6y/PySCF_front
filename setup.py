"""
PySCF_Front セットアップスクリプト
"""
from setuptools import setup, find_packages
from pathlib import Path

# README の読み込み
readme_path = Path(__file__).parent / "README.md"
long_description = readme_path.read_text(encoding="utf-8") if readme_path.exists() else ""

setup(
    name="pyscf-front",
    version="0.0.1",
    description="PySide6-based GUI for PySCF quantum chemistry calculations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="PySCF_Front Development Team",
    author_email="dev@pyscf-front.org",
    url="https://github.com/your-org/pyscf-front",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.13",
        "Programming Language :: Python :: 3.14",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    python_requires=">=3.13",
    install_requires=[
        "PySide6>=6.7.0",
        "PySide6-Addons>=6.7.0",
        "PySide6-Essentials>=6.7.0",
        "pyscf>=2.9.0",
        "numpy>=1.26.0",
        "scipy>=1.13.0",
        "matplotlib>=3.9.0",
        "vtk>=9.3.0",
        "sqlalchemy>=2.0.30",
        "loguru>=0.7.0",
        "pyyaml>=6.0.0",
        "python-dotenv>=1.0.0",
    ],
    extras_require={
        "dev": [
            "pytest>=8.2.0",
            "pytest-qt>=4.4.0",
            "pytest-cov>=5.0.0",
            "black>=24.4.0",
            "flake8>=7.0.0",
            "mypy>=1.10.0",
        ],
        "gpu": [
            "gpu4pyscf-cuda12x>=1.4.0",
        ],
        "mcp": [
            "mcp>=1.0.0",
            "fastmcp>=0.14.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "pyscf-front=pyscf_front.main:main",
            "pyscf-front-dev=run_dev:main",
        ],
    },
    package_data={
        "pyscf_front": [
            "resources/icons/*",
            "resources/styles/*",
            "resources/ui_files/*",
            "translations/*.qm",
        ],
    },
    include_package_data=True,
    zip_safe=False,
)
from pathlib import Path

from setuptools import setup

here = Path(__file__).resolve().parent


README = (here / "README.md").read_text(encoding="utf-8")

requirements = [
    line.strip()
    for line in (here / "requirements.txt").read_text(encoding="utf-8").split("\n")
    if line.strip()
]


setup(
    name="fuchsian-toolkit",
    version="1.0",
    packages=["fuchsian_toolkit"],
    description="Tools to constract rational transformation matrix",
    include_package_data=True,
    long_description=README,
    long_description_content_type="text/markdown",
    author="Giorgi Kakulashvili",
    url="https://github.com/giorgi94/sc-tools",
    keywords=["fuchsian", "rational matrix"],
    platforms=["OS Independent"],
    install_requires=requirements,
)

from setuptools import setup, find_packages


with open("README.md", "r", encoding="utf-8") as fin:
    long_description = fin.read()

setup(
    name="navann-pkg-mapqzero",
    version="0.0.1",
    author="Andrei Rajkovic",
    author_email="rajkovic.1@osu.edu",
    description="A small variant annotations package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/andreirajkovic/navann",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=find_packages(),
    python_requires='>=3.6',
    entry_points={
        'console_scripts': [
            'NAVANN = VariantFileAnnotator:main',
        ],
    },    
)
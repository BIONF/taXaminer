from setuptools import setup, find_packages

with open("README.md", "r") as input:
    long_description = input.read()

setup(
    name="taxaminer",
    version="0.7.0",
    python_requires='>=3.7.0',
    description="Interactive exploration of biodiverse genome assemblies",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Freya Arthen",
    author_email="f.arthen@bio.uni-frankfurt.de",
    url="https://github.com/BIONF/taXaminer",
    packages=find_packages(),
    package_data={'': ['*']},
    install_requires=[
        'beautifulsoup4',
        'jsmin',
        'PyYAML',
        'biopython',
        'numpy',
        'pandas',
        'plotly',
        'requests',
        'taxopy',
        'scikit-learn',
        'scipy',
        'umap-learn',
        'kaleido',
        'pysam'
    ],
    entry_points={
        'console_scripts': ["taxaminer.run = taxaminer.runTaXaminer:main",
                            "taxaminer.setup = taxaminer.setupTaXaminer:main",
                            ],
    },
    license="MIT",
    classifiers=[
        "Environment :: Console",
        "Intended Audience :: End Users/Desktop",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
    ],
)

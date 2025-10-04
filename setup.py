from setuptools import setup, find_packages

setup(
    name='ecdnainspector',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'pyyaml',
        'pandas',
        'numpy==1.24',
        'seaborn',
        'matplotlib',
        'scikit-learn'
    ],
    entry_points={
        'console_scripts': [
            'ecI=ecDNAInspector.cli:main',
        ],
    },
    author='Sophia J. Pribus',
    description='ecDNAInspector: A tool for filtering and analyzing extrachromosomal DNA',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7',
)
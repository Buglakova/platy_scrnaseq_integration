from setuptools import setup, find_packages

setup(
    name='cluster_map',
    version='1.0',
    description='Package for mapping clusters annotated in scRNAseq to spatial transcriptomics data',
    url='https://github.com/Buglakova/platy_scrnaseq_integration',
    packages=find_packages(include=['cluster_map']),
    python_requires='>=3.6',
    install_requires=[],
    author='Elena Buglakova, Cyril Cros',
    author_email='elena.buglakova@embl.de',
    license='MIT'
)

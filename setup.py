from setuptools import setup, find_packages, find_namespace_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
        name='mirge',
        version='0.1.46',
        author='Arun Patil and Marc Halushka',
        author_email='mhalush1@jhmi.edu',
        url='https://test.pypi.org/legacy/',
        description='Comprehensive analysis of small RNA sequencing Data', 
        long_description=long_description,
        keywords=['miRge', 'small RNA analysis', 'NGS', 'bioinformatics tools'],  # arbitrary keywords
        license='MIT',
        package_dir={'mirge': 'mirge'},
        packages=find_packages(),
        package_data = {'':['models/*.pkl', 'models/*.txt', 'rScripts/*.R']},
        install_requires=['biopython==1.77','cutadapt==2.7','future==0.18.2','joblib==0.15.1','matplotlib==3.2.1',
            'numpy==1.18.4','pandas==1.0.3','reportlab==3.5.42','scikit-learn==0.23.1','scipy==1.4.1',
            'umi-tools==1.0.1'],
        entry_points={'console_scripts': ['miRge3.0 = mirge.__main__:main']},
        classifiers=[
            "Development Status :: 1 - Planning",
            "Environment :: Console",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Natural Language :: English",
            "Programming Language :: Python :: 3.8",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
            ],
        include_package_data=True,
        python_requires='>=3.8',
)

import setuptools

setuptools.setup(
    name="singlecell",
    packages=setuptools.find_packages(),
    python_requires='>=3.6',
    scripts=[
        'mt_prepend_gtf.py',
        'create_cellranger_ref.py',
    ],
    package_data={
        'mkref': ['create_cellranger_ref.sh'],
    },
)

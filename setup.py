from setuptools import setup, find_packages, Extension
setup(name='Mars',      
version='1.2',      
description='Multiple Alignment-based Refinement of SVs',      
author='Xin Zhou',      
author_email='maizie.zhou@vanderbilt.edu',      
install_requires=['pysam','pandas','minimap2','openpyxl','samtools'],         
license='MIT',      
packages=find_packages(),
package_data={'bin' : ['paftools/*',
                     'k8-0.2.4/*',
                     'trf_tools/*',
                     ]},  
entry_points={'console_scripts':['MARS_step1=bin.MARS_step1:main',
'MARS_step2=bin.MARS_step2:main'
               ]},
zip_safe=False)

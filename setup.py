from setuptools import setup
setup(name='Mars',      
version='1.0',      
description='Multiple Alignment-based Refinement of SVs',      
author='Xin Zhou',      
author_email='mazie.zhou@vanderbilt.edu',      
install_requires=['pysam','pandas','minimap2','openpyxl','samtools'],         
license='MIT',      
packages=['bin',],
package_data={'bin' : ['paftools/*',
                     'k8-0.2.4/*',
                     'trf_tools/*',
                     ]},  
entry_points={'console_scripts':['MARS_step1=bin.MARS_step1:main',
'MARS_step2=bin.MARS_step2:main'
               ]},
zip_safe=False)

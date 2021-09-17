# change permission
cd bin
chmod +x *.py
chmod +x muscle3.8.31_i86linux64
chmod -R 777 k8-0.2.4
cd ..

if ! [ -x "$(command -v samtools)" ];
then
    echo 'Error: samtools is not installed...'
    echo 'using Conda to install...'
    conda install -c bioconda samtools
else
    echo 'using existing samtools...'
fi

if ! [ -x "$(command -v minimap2)" ];
then
    echo 'Error: minimap2 is not installed...'
    echo 'using Conda to install...'
    conda install -c bioconda minimap2
else
    echo 'using existing minimap2...'
fi



echo 'You have installed Aquila dependencies and downloaded the source files successfully!'
 




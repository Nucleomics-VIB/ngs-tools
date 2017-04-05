[(Nucleomics-VIB)](https://github.com/Nucleomics-VIB)
![ngs-tools](pictures/pacbio_icon.png) - NGS-Tools
==========

*All tools presented below have only been tested by me and may contain bugs, please let me know if you find some. Each tool relies on dependencies normally listed at the top of the code (cpan for perl and cran for R will help you add them)*

Please refer to the accompanying **[wiki](https://github.com/Nucleomics-VIB/ngs-tools/wiki)** for examples and workflows.

### **makeENSref.sh**

The bash file **[makeENSref.sh](/ngs-tools/makeENSref.sh)** creates a series of reference files from ENSembl FTP downloads.
```bash
# Usage: makeENSref.sh
# -o <organism (default to <homo_sapiens>)> 
# -b <build number (default to <GRCh38>)> 
# -r <release number (default to 88)>
# script version 1.0, 2017_04_05
# [-h for this help]
```
<hr>

<h4>Please send comments and feedback to <a href="mailto:nucleomics.bioinformatics@vib.be">nucleomics.bioinformatics@vib.be</a></h4>

<hr>

![Creative Commons License](http://i.creativecommons.org/l/by-sa/3.0/88x31.png?raw=true)

This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/).

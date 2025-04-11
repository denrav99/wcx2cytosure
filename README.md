# wcx2cytosure

This tool converts wisecondorx output to the “.CGH” format used by the
commercial
[CytoSure Interpret Software](https://www.ogt.com/products/246_cytosure_interpret_software)
by OGT (Oxford Gene Technology). CytoSure is made for displaying oligo array measurements.
It works on a set of probes, which this tool emulates. 

## Usage
Requires an input .bins and .aberrations file from wisecondorx. 

    INSTALL:
        git  clone https://github.com/denrav99/wcx2cytosure.git
        cd wcx2cytosure
        pip install -e .
	
    or use Singularity  
    
    	singularity pull --arch amd64 library://ravinale/wcx2cytosure/wcx2cytosure2:latest

    RUNNING:
        
    wcx2cytosure --wisecondorx_cov  <input.bins.bed> --wisecondorx_aberrations <input.aberrations.bed> --out <output.cgh> --wcx_size <smallest aberration size (int)> (optional)
    
## Notes on the file format

- CGH is in XML format. The company does not seem to have a schema file.
- The parser in CytoSure accepts reformatted XML files. This is ok:

      xmllint --format file.cgh > pretty.xml
- One probe can be made up of more than one spot:

      <probe ...>
        <spot ... />
        <spot ... />
      </probe>

- The probes in the file do not have to be sorted by chromosome and/or
  coordinate.
- There are `<probe>` elements without sequence, without chromosome name:

      <probe name="probename" sequence="">
        <spot index="56774" row="334" column="164" red="123.4" green="345.6" gSNR="1.0" rSNR="50.0" outlier="false"/>
      </probe>

- Coordinates are 1-based (since stop - start = 59 and length of sequence is 60)
- Removing the 'sequence' attribute does not work, but setting it to a fake one does
- The log2 ratio is computed as log2(green/red)
- Normalization segments look like this:

      <segment chrId="1" numProbes="10" start="7000" stop="1000000" average="-0.003"/>
- Segments can overlap each other
- Names "X" and "Y" are not used for the corresponding chromosomes. Instead,
  the `X` chromosome is stored as `23`, they `Y` chromosome as `24` (in probes,
  segments and aberrations).


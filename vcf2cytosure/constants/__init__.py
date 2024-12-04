# constants

PROBE_SPACING = 100000
MAX_HEIGHT = 4
MIN_HEIGHT = -4
ABERRATION_HEIGHTS = {
	'DEL': -1.0,
	'DUP': +1.0,
}

CONTIG_LENGTHS_37 = {
	'1':  249250621,
	'2':  243199373,
	'3':  198022430,
	'4':  191154276,
	'5':  180915260,
	'6':  171115067,
	'7':  159138663,
	'8':  146364022,
	'9':  141213431,
	'10': 135534747,
	'11': 135006516,
	'12': 133851895,
	'13': 115169878,
	'14': 107349540,
	'15': 102531392,
	'16':  90354753,
	'17':  81195210,
	'18':  78077248,
	'19':  59128983,
	'20':  63025520,
	'21':  48129895,
	'22':  51304566,
	'X':  155270560,
	'Y':   59373566,
}

CONTIG_LENGTHS_38 = {
	'1':    248956422,
	'2':    242193529,
	'3':	198295559,
	'4':    190214555,
	'5':    181538259,
	'6':    170805979,
	'7':    159345973,
	'8':    145138636,
	'9':    138394717,
	'10':   133797422,
	'11':   135086622,
	'12':   133275309,
	'13':   114364328,
	'14':   107043718,
	'15':   101991189,
	'16':   90338345,
	'17':   83257441,
	'18':   80373285,
	'19':   58617616,
	'20':   64444167,
	'21':   46709983,
	'22':   50818468,
	'X':    15604089,
	'Y':    57227415,
}

CHROM_RENAME = {'X': '23', 'Y': '24'}

CGH_TEMPLATE_37 = u"""
<data formatVersion="2">
<pgdData><pgdDataEntry key="SPECIMEN_TYPE" value="BLASTOMERE"/></pgdData>
<noResults>
</noResults>
<pgd_reagents/>
<cgh mother="-1" father="-1" genomeBuild="hg19" softwareVersion="4.8.32" batched="false">
  <submission design="031035" feFile="{}.txt" cghFile="{}.cgh" scanDate="1462520414000" barcode="{}" sampleCy3="true">
  <notes/>
  <sample sampleId="{}" male="{}"><phenotype/></sample>
  <reference sampleId="Promega {}" male="{}"><phenotype/></reference>
  <extra>
    <datum category="Nanodrop" type="Sample DNA (ng)" dataType="Float"/>
    <datum category="Sample Extraction" type="Sample Arrival Date" dataType="Date"/>
    <datum category="Sample Extraction" type="Specimen Type" dataType="List">Blood</datum>
    <datum category="Labelling" type="Lab User" dataType="String"/>
    <datum category="Hyb &amp; Wash" type="Cot1 Batch" dataType="String"/>
    <datum category="Sample Extraction" type="Extraction Method" dataType="String"/>
    <datum category="General" type="Reference Concentration" dataType="Float"/>
    <datum category="Sample Extraction" type="A260/A280" dataType="Float"/>
    <datum category="General" type="Assigned Technologist" dataType="String">MF</datum>
    <datum category="Nanodrop" type="Reference DNA (pmoles)" dataType="Float"/>
    <datum category="Labelling" type="Columns Used" dataType="List">Qiagen</datum>
    <datum category="Nanodrop" type="Sample DNA (A260/A280)" dataType="Float"/>
    <datum category="Labelling" type="Column Batch No" dataType="String"/>
    <datum category="Sample Extraction" type="Extracted By" dataType="String"/>
    <datum category="Hyb &amp; Wash" type="Hyb Protocol" dataType="String"/>
    <datum category="Labelling" type="Experiment Date" dataType="Date"/>
    <datum category="General" type="Case Status" dataType="List"/>
    <datum category="Labelling" type="Lab Notes" dataType="Text"/>
    <datum category="Hyb &amp; Wash" type="Hyb Buffer Batch" dataType="String"/>
    <datum category="Nanodrop" type="Reference DNA (A260/A280)" dataType="Float"/>
    <datum category="Labelling" type="Labelling Reagent Batch No" dataType="String"/>
    <datum category="Hyb &amp; Wash" type="Wash Protocol" dataType="String"/>
    <datum category="Labelling" type="Labelling Protocol" dataType="List">Enzo</datum>
    <datum category="Nanodrop" type="Reference DNA (ng)" dataType="Float"/>
    <datum category="Nanodrop" type="Sample DNA (pmoles)" dataType="Float"/>
  </extra>
</submission>
<excludedRegions>
</excludedRegions>
<qc>
<aqf key="SPIKES" value="null"/>
<aqf key="DLRSPREAD" value="0.1"/>
<aqf key="RED_SIGNAL_INTENSITY" value="1000.0"/>
<aqf key="GREEN_SIGNAL_INTENSITY" value="1000.0"/>
<aqf key="BG_NOISE_RED" value="5.0"/>
<aqf key="BG_NOISE_GREEN" value="5.0"/>
<aqf key="RED_SNR" value="100.0"/>
<aqf key="GREEN_SNR" value="100.0"/>
<aqf key="COMPARATIVE_SIGNAL_INTENSITY" value="0.11111111111"/>
<aqf key="REPRODUCIBILITY_GREEN" value="0.1111111"/>
<aqf key="REPRODUCIBILITY_RED" value="0.1111111"/>
<aqf key="NEG_CONTROL_RED" value="1.1111111"/>
<aqf key="NEG_CONTROL_GREEN" value="1.1111111"/>
<aqf key="FLAG_PERC" value="0.00722363"/>
<aqf key="SAT_PERC" value="0.0"/>
<aqf key="WAVINESS" value="0.0011111"/>
<aqf key="SNP_TROUGH_PEAK_RATIO" value="null"/>
<aqf key="SNP_GREY_AREA" value="null"/>
<aqf key="SNP_Red_Signal_Intensity" value="NaN"/>
<aqf key="SNP_Green_Signal_Intensity" value="NaN"/>
<aqf key="SNP_Ratio_Separation" value="0.0"/>
<aqf key="Percentage_Homozygosity" value="NaN"/>
<aqf key="SD" value="0.111111"/>
</qc>
<probes></probes>
<segmentation type="NORMALIZED"></segmentation>
</cgh>
</data>
"""

CGH_TEMPLATE_38 = u"""
<data formatVersion="2">
<pgdData><pgdDataEntry key="SPECIMEN_TYPE" value="BLASTOMERE"/></pgdData>
<noResults>
</noResults>
<pgd_reagents/>
<cgh mother="-1" father="-1" genomeBuild="hg38" softwareVersion="4.10.41" batched="false">
  <submission design="031035" feFile="{}.txt" cghFile="{}.cgh" scanDate="1462520414000" barcode="{}" sampleCy3="true">
  <notes/>
  <sample sampleId="{}" male="{}"><phenotype/></sample>
  <reference sampleId="Promega {}" male="{}"><phenotype/></reference>
  <extra>
    <datum category="Nanodrop" type="Sample DNA (ng)" dataType="Float"/>
    <datum category="Sample Extraction" type="Sample Arrival Date" dataType="Date"/>
    <datum category="Sample Extraction" type="Specimen Type" dataType="List">Blood</datum>
    <datum category="Labelling" type="Lab User" dataType="String"/>
    <datum category="Hyb &amp; Wash" type="Cot1 Batch" dataType="String"/>
    <datum category="Sample Extraction" type="Extraction Method" dataType="String"/>
    <datum category="General" type="Reference Concentration" dataType="Float"/>
    <datum category="Sample Extraction" type="A260/A280" dataType="Float"/>
    <datum category="General" type="Assigned Technologist" dataType="String">MF</datum>
    <datum category="Nanodrop" type="Reference DNA (pmoles)" dataType="Float"/>
    <datum category="Labelling" type="Columns Used" dataType="List">Qiagen</datum>
    <datum category="Nanodrop" type="Sample DNA (A260/A280)" dataType="Float"/>
    <datum category="Labelling" type="Column Batch No" dataType="String"/>
    <datum category="Sample Extraction" type="Extracted By" dataType="String"/>
    <datum category="Hyb &amp; Wash" type="Hyb Protocol" dataType="String"/>
    <datum category="Labelling" type="Experiment Date" dataType="Date"/>
    <datum category="General" type="Case Status" dataType="List"/>
    <datum category="Labelling" type="Lab Notes" dataType="Text"/>
    <datum category="Hyb &amp; Wash" type="Hyb Buffer Batch" dataType="String"/>
    <datum category="Nanodrop" type="Reference DNA (A260/A280)" dataType="Float"/>
    <datum category="Labelling" type="Labelling Reagent Batch No" dataType="String"/>
    <datum category="Hyb &amp; Wash" type="Wash Protocol" dataType="String"/>
    <datum category="Labelling" type="Labelling Protocol" dataType="List">Enzo</datum>
    <datum category="Nanodrop" type="Reference DNA (ng)" dataType="Float"/>
    <datum category="Nanodrop" type="Sample DNA (pmoles)" dataType="Float"/>
  </extra>
</submission>
<excludedRegions>
</excludedRegions>
<qc>
<aqf key="SPIKES" value="null"/>
<aqf key="DLRSPREAD" value="0.1"/>
<aqf key="RED_SIGNAL_INTENSITY" value="1000.0"/>
<aqf key="GREEN_SIGNAL_INTENSITY" value="1000.0"/>
<aqf key="BG_NOISE_RED" value="5.0"/>
<aqf key="BG_NOISE_GREEN" value="5.0"/>
<aqf key="RED_SNR" value="100.0"/>
<aqf key="GREEN_SNR" value="100.0"/>
<aqf key="COMPARATIVE_SIGNAL_INTENSITY" value="0.11111111111"/>
<aqf key="REPRODUCIBILITY_GREEN" value="0.1111111"/>
<aqf key="REPRODUCIBILITY_RED" value="0.1111111"/>
<aqf key="NEG_CONTROL_RED" value="1.1111111"/>
<aqf key="NEG_CONTROL_GREEN" value="1.1111111"/>
<aqf key="FLAG_PERC" value="0.00722363"/>
<aqf key="SAT_PERC" value="0.0"/>
<aqf key="WAVINESS" value="0.0011111"/>
<aqf key="SNP_TROUGH_PEAK_RATIO" value="null"/>
<aqf key="SNP_GREY_AREA" value="null"/>
<aqf key="SNP_Red_Signal_Intensity" value="NaN"/>
<aqf key="SNP_Green_Signal_Intensity" value="NaN"/>
<aqf key="SNP_Ratio_Separation" value="0.0"/>
<aqf key="Percentage_Homozygosity" value="NaN"/>
<aqf key="SD" value="0.111111"/>
</qc>
<probes></probes>
<segmentation type="NORMALIZED"></segmentation>
</cgh>
</data>
"""

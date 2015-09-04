# panel is a dict containing all of the data structures that define a panel
# 	Amplicons:		genomic range
# 	Exons:			genomic range
# 	Transcripts:		genomic range
# 	Depth:
#	AllTranscripts: 	genomic range
#	AllExons:		genomic range
#	Excluded:		list of excluded amplicons
#	Variants_Disease:	genomic range
#	Variants_Gene:		genomic range
#	Variants_Mutation:	genomic range
#	Filenames		dict of all files in the panel directory with filetype as the key
#


import os, re, pdb
import gr


class CoverMiException(Exception):
    pass

filetypes = ["Excluded", "Reference", "Depth", "Targets", "Targets", "Manifest", "Variants", "Canonical", "DesignStudio", "Disease_Names", "SNPs"]
regexps = ["^chr[1-9XYM][0-9]?:[0-9]+-[0-9]+\\s*\n", #Single column of amplicon names in the format chr1:12345-67890 - excluded
           "^.+?\t.+?\tchr.+?\t[+-]\t[0-9]+\t[0-9]+\t[0-9]+\t[0-9]+\t[0-9]+\t[0-9,]+\t[0-9,]+\n", #Eleven columns - refflat
           "^[0-9]+\\s*\n", #Single number - depth
           "^[a-zA-Z][a-zA-Z0-9-]*\\s*\n", #Single column - targets
           "^[a-zA-Z][a-zA-Z0-9-]*\\s+[a-zA-Z0-9_]+\\s*\n", #Single column - targets
           "^\\[Header\\]", #Manifest
           "^.+?\t.+?\t.+?\t[a-zA-Z0-9]+\t(null|chr[1-9XYM][0-9]?)\t(null|[0-9]+)\t(null|[0-9]+)\t(null|[+-])\t", #Nine columns - variants
           "^[a-zA-Z][a-zA-Z0-9]*\t[a-zA-Z][a-zA-Z0-9-]*\\s+[a-zA-Z0-9_]+\\s*\n", #Two columns - canonical
           "^chr[0-9XYM]+\t[0-9]+\t[0-9]+\t[^\t]+\t[0-9]+\t[+-]\n", #Six column bedfile - design
           "^#Variants Disease Name Translation\n", #Disease_Names
           "^chr[1-9XYM][0-9]?:[0-9]+ [ATCG\.]+>[ATCG\.]\n"] #SNPs

def identify(panel_path):
    path = {}
    for root, dirnames, filenames in os.walk(panel_path):
        for filename in filenames:
            filepath = os.path.join(root, filename)
            with file(filepath, "rU") as f:
                testlines = [f.readline(), f.readline()]
            
            alreadyfound = ""
            for testline in testlines:
                if not testline.endswith("\n"):
                    testline = testline+"\n"
                for regexp, filetype in zip(regexps, filetypes):
                    if re.match(regexp, testline) is not None:
                        if alreadyfound != "":
                            raise CoverMiException("Unable to uniquely identify file {0}, matches both {1} and {2} file formats".format(filepath, filetype, alreadyfound))
                        if filetype in path:
                            raise CoverMiException("Files {0} and {1} both match {2} file format".format(filepath, path[filetype], filetype))
                        alreadyfound = filetype
                        path[filetype] = filepath
                if alreadyfound != "":
                    break
    return path


def load(panel_path, design=False):
    panel = {}    
    path = identify(panel_path)

    if "Reference" not in path:
       raise CoverMiException("Unable to identify a Reference file in {0}".format(panel_path))
    if design and "Manifest" not in path and "DesignStudio" not in path:
        raise CoverMiException("Unable to identify a Manifest file or a Design Studio Bedfile in {0}".format(panel_path))
    if "Manifest" in path and "DesignStudio" in path:
        raise CoverMiException("Both Manifest and Design Studio Files found in {0}".format(panel_path))

    if "Depth" in path:
        with file(path["Depth"], "rU") as f:
            depth = f.read().rstrip()
        if not depth.isdigit():
            raise CoverMiException("Depth file is not numeric")
        panel["Depth"] = int(depth)    

    if "Manifest" in path:
        excluded = []
        print "Loading manifest file: {0}".format(os.path.basename(path["Manifest"]))
        if "Excluded" in path:
            with file(path["Excluded"], "rU") as f:
                for line in f:
                    line = line.rstrip()
                    if line != "":
                        excluded.append(line)
            panel["Excluded"] = excluded
            if design:
                panel["AllAmplicons"] = gr.load_manifest(path["Manifest"], [])
        panel["Amplicons"] = gr.load_manifest(path["Manifest"], panel["Excluded"] if ("Excluded" in panel) else [])
    elif "DesignStudio" in path:
        print "Loading Design Studio amplicons bedfile: {0}".format(os.path.basename(path["DesignStudio"]))
        panel["Amplicons"] = gr.load(path["DesignStudio"])
    elif "Targets" in path:
        print "No manifest or design studio bedfile in panel. Will perform coverage analysis over exons of genes in targets file"
    else:
        raise CoverMiException("WARNING. No manifest or design studio bedfile or targeted genes file in panel. Not yet set up to perform whole exome coverage analysis.")

    if "Reference" in path:
        print "Loading reference file: {0}".format(os.path.basename(path["Reference"]))
        if "Targets" in path:
            print "Loading targeted genes file: {0}".format(os.path.basename(path["Targets"]))
        if "Canonical" in path:
            print "Loading canonical gene list: {0}".format(os.path.basename(path["Canonical"]))
        else:
            print "WARNING. No canonical gene list found. Loading all transcripts of targeted genes"
        panel["Exons"], panel["Transcripts"] = gr.load_refflat(path["Reference"], 
                                                                         path["Targets"] if ("Targets" in path) else "", 
                                                                         path["Canonical"] if ("Canonical" in path) else "")

        if "Targets" not in path and "Amplicons" in panel:
            print "WARNING. No file found identifying genes of interest. All genes touched by an amplicon will be assumed to be of interest"
            panel["Transcripts"] = panel["Transcripts"].touched_by(panel["Amplicons"])
            panel["Exons"] = panel["Exons"].touched_by(panel["Transcripts"]).subset2(panel["Transcripts"].names)
        if design:
            panel["AllExons"], panel["AllTranscripts"] = gr.load_refflat(path["Reference"], "", path["Canonical"] if ("Canonical" in path) else "")

    if "Variants" in path:
        print "Loading variants file: {0}".format(os.path.basename(path["Variants"]))
        if "Disease_Names" in path:
            print "Loading disease names file: {0}".format(os.path.basename(path["Disease_Names"]))
            disease_names_path = path["Disease_Names"]
        else:
            disease_names_path = None
        panel["Variants_Disease"] = gr.load_variants(path["Variants"], "disease", disease_names=disease_names_path)
        panel["Variants_Gene"] = gr.load_variants(path["Variants"], "gene", genes_of_interest=panel["Transcripts"], disease_names=disease_names_path)
        panel["Variants_Mutation"] = gr.load_variants(path["Variants"], "mutation", genes_of_interest=panel["Transcripts"], disease_names=disease_names_path)

    panel["Filenames"] = { "Panel" : os.path.basename(panel_path.rstrip(os.pathsep)) }
    for filetype in path:
        panel["Filenames"][filetype] = os.path.splitext(os.path.basename(path[filetype]))[0]
        
    return panel


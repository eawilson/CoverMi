# panel is a dict containing all of the data structures that define a panel
# 	Amplicons:		genomic range
# 	Exons:			genomic range
# 	Transcripts:		genomic range
# 	Depth:
#	Variants_Disease:	genomic range
#	Variants_Gene:		genomic range
#	Variants_Mutation:	genomic range
#	Filenames:		dict of all files in the panel directory with filetype as the key
#       Options:		dict of all options, including depth

# the following are only included for a design panel
#	 AllTranscripts: 	genomic range
#	 AllExons:		genomic range
#	 Excluded:		list of excluded amplicons


import os, re, pdb
from gr import Gr


class CoverMiException(Exception):
    pass


filetypes = ["Excluded", "Reference", "Depth", "Targets", "Targets", "Manifest", "Variants", 
             "Canonical", "Canonical", "DesignStudio", "Disease_Names", "SNPs", "VCF", "Options"]
regexps = ["^chr[1-9XYM][0-9]?:[0-9]+-[0-9]+\\s*\n", #Single column of amplicon names in the format chr1:12345-67890 - excluded
           "^.+?\t.+?\tchr.+?\t[+-]\t[0-9]+\t[0-9]+\t[0-9]+\t[0-9]+\t[0-9]+\t[0-9,]+\t[0-9,]+\\s*\n", #Eleven columns - refflat
           "^[0-9]+\\s*\n", #Single number - depth
           "^[a-zA-Z][a-zA-Z0-9-]*\\s*\n", #Single column - targets
           "^[a-zA-Z][a-zA-Z0-9-]* +[a-zA-Z0-9_]+\\s*\n", #Single column - targets
           "^\\[Header\\]", #Manifest
           "^.+?\t.+?\t.+?\t[a-zA-Z0-9]+\t(null|chr[1-9XYM][0-9]?)\t(null|[0-9]+)\t(null|[0-9]+)\t(null|[+-])\t", #Nine columns - variants
           "^[a-zA-Z][a-zA-Z0-9]*\t[a-zA-Z][a-zA-Z0-9-]* +[a-zA-Z0-9_]+\\s*\n", #Two columns - canonical
           "^[a-zA-Z][a-zA-Z0-9]*\t[a-zA-Z0-9_]+\\s*\n", #Two columns - canonical
           "^chr[0-9XYM]+\t[0-9]+\t[0-9]+\t[^\t]+\t[0-9]+\t[+-]\\s*\n", #Six column bedfile - design
           "^#Variants Disease Name Translation\\s*\n", #Disease_Names
           "^chr[1-9XYM][0-9]?:[0-9]+ [ATCG\.]+>[ATCG\.]\\s*\n", #SNPs
           "^##fileformat=VCFv[0-9\.]+\\s*\n", #VCF
           "#CoverMi options\\s*\n"] # CoverMi panel options


class Panel(dict):


    def __init__(self, panel_path):
        self.panel_path = panel_path
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
                                raise CoverMiException("Not a valid CoverMi panel\nUnable to uniquely identify file {0}, matches both {1} and {2} file formats".format(filepath, filetype, alreadyfound))
                            if filetype in self:
                                raise CoverMiException("Not a valid CoverMi panel\nFiles {0} and {1} both match {2} file format".format(filepath, self[filetype], filetype))
                            alreadyfound = filetype
                            self[filetype] = filepath
                    if alreadyfound != "":
                        break
        return


    def read_options(self):
        options = {}
        if "Depth" in self:
            with file(self["Depth"], "rU") as f:
                options["Depth"] = f.read().strip()
        if "Options" in self:
            print "Loading options file: {0}".format(os.path.basename(self["Options"]))
            with file(self["Options"], "rU") as f:
                for line in f:
                    line = line.strip()
                    if line != "" and not line.startswith("#"):
                        line = line.split("=")
                        if len(line) != 2:
                            raise CoverMiException("Malformed options file: {0}".format("=".join(line)))
                        if line[0] in options:
                            raise CoverMiException("Duplicate option: {0}".format(line[0]))
                        options[line[0]] = line[1]
        if "Depth" in options:
            if  options["Depth"].isdigit():
                options["Depth"] = int(options["Depth"])
            else:
                raise CoverMiException("Depth is not numeric: {0}".format(panel["Options"]["Depth"]))
        return options


    def write_options(self, options):
        if "Options" not in self:
            self["Options"] = os.path.join(self.panel_path, "covermi_options.txt")
            if os.path.lexists(self["Options"]):
                raise CoverMiException("File covermi_options.txt exists but is not of correct format")
        with file(self["Options"], "wt") as f: 
            f.write("#CoverMi options\n")
            for key, value in options.iteritems():
                f.write("{0}={1}\n".format(key, value))
        if "Depth" in self:
            os.unlink(self["Depth"])


    def load(self, design=False):
        panel = {}    

        if "Reference" not in self:
            raise CoverMiException("Unable to identify a Reference file in {0}".format(self.panel_path))
        if design and "Manifest" not in self and "DesignStudio" not in self:
            raise CoverMiException("Unable to identify a Manifest file or a Design Studio Bedfile in {0}".format(self.panel_path))
        if "Manifest" in self and "DesignStudio" in self:
            raise CoverMiException("Both Manifest and Design Studio Files found in {0}".format(self.panel_path))

        panel["Options"] = self.read_options()

        if "Manifest" in self:
            print "Loading manifest file: {0}".format(os.path.basename(self["Manifest"]))
            panel["Amplicons"] = Gr.load_manifest(self["Manifest"])
            if "Excluded" in self:
                excluded = []
                with file(self["Excluded"], "rU") as f:
                    for line in f:
                        line = line.strip()
                        if line != "" and not line.startswith("#"):
                            excluded.append(line)
                if design:
                    panel["ExcludedAmplicons"] = panel("Amplicons").subset(excluded)
                panel["Amplicons"] = panel("Amplicons").subset(excluded, exclude=True)

        elif "DesignStudio" in self:
            print "Loading Design Studio amplicons bedfile: {0}".format(os.path.basename(self["DesignStudio"]))
            panel["Amplicons"] = Gr.load_bedfile(self["DesignStudio"])
        elif "Targets" in self:
            print "No manifest or design studio bedfile in panel. Will perform coverage analysis over all exons of genes in targets file"
        else:
            raise CoverMiException("WARNING. No manifest or design studio bedfile or targeted genes file in panel. Aborting.")

        if "Reference" in self:
            print "Loading reference file: {0}".format(os.path.basename(self["Reference"]))
            if "Targets" in self:
                print "Loading targeted genes file: {0}".format(os.path.basename(self["Targets"]))
            if "Canonical" in self:
                print "Loading canonical gene list: {0}".format(os.path.basename(self["Canonical"]))
            else:
                print "WARNING. No canonical gene list found. Loading all transcripts of targeted genes"
            panel["Exons"], panel["Transcripts"] = Gr.load_refflat(self["Reference"], 
                                                                   file(self["Targets"], "rU") if ("Targets" in self) else None, 
                                                                   file(self["Canonical"], "rU") if ("Canonical" in self) else None)

            if "Targets" not in self and "Amplicons" in panel:
                print "WARNING. No file identifying targeted genes. All genes touched by an amplicon will be assumed to be of interest"
                panel["Transcripts"] = panel["Transcripts"].touched_by(panel["Amplicons"])
                panel["Exons"] = panel["Exons"].touched_by(panel["Transcripts"]).subset(panel["Transcripts"].names)
            if design:
                panel["AllExons"], panel["AllTranscripts"] = Gr.load_refflat(self["Reference"], None, file(self["Canonical"], "rU") if ("Canonical" in self) else "")

        if "Variants" in self:
            print "Loading variants file: {0}".format(os.path.basename(self["Variants"]))
            if "Disease_Names" in self:
                print "Loading disease names file: {0}".format(os.path.basename(self["Disease_Names"]))
            panel["Variants_Disease"] = Gr.load_variants(self["Variants"], "disease",
                                                         disease_names=file(self["Disease_Names"], "rU") if ("Disease_Names" in self) else None)
            panel["Variants_Gene"] = Gr.load_variants(self["Variants"], "gene", genes_of_interest=panel["Transcripts"].names, 
                                                      disease_names=file(self["Disease_Names"], "rU") if ("Disease_Names" in self) else None)
            panel["Variants_Mutation"] = Gr.load_variants(self["Variants"], "mutation", genes_of_interest=panel["Transcripts"].names,
                                                          disease_names=file(self["Disease_Names"], "rU") if ("Disease_Names" in self) else None)
            
        panel["Filenames"] = { "Panel" : os.path.basename(self.panel_path.rstrip(os.pathsep)) }
        for filetype in self:
            panel["Filenames"][filetype] = os.path.splitext(os.path.basename(self[filetype]))[0]
        
        return panel


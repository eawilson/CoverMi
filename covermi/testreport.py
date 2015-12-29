import pkg_resources
try:
    from covermiglobals import *
    import GenomicRange 
except ImportError:
    from gr import Gr as GenomicRange
    VERSION = pkg_resources.require("CoverMi")[0].version

import time, os, pdb


def create(coverage, panel, output_path):

    try:
        output_file = os.path.join(output_path, panel["Filenames"]["Output"]+"_covermi_test_report.txt")
    except KeyError:
        output_file = output_path+"_covermi_test_report.txt"
    with file(output_file, "wt") as f:
        f.write("{0}{1}\n".format("CoverMi version:", VERSION))
        
        try:
            minimum_depth = panel["Depth"]
        except KeyError:
            minimum_depth = panel["Options"]["Depth"]

        targeted_exons = panel["Exons"] if ("Exons" in panel) else panel["Transcripts"]
        if "Amplicons" in panel:
            targeted_range = panel["Amplicons"].merged
            targeted_exons = targeted_exons.overlapped_by(targeted_range)
        else:
            targeted_range = panel["Transcripts"]
            targeted_exons = targeted_exons

        try:
            i = coverage.info(targeted_exons.renamed("all"), minimum_depth)[0]
        except AttributeError:
            i = coverage.calculate(targeted_exons, minimum_depth, total=True)
        f.write("\n{0} of {1} bases ({2:0.0f}%) covered at depth >= {3}  (mean depth {4})\n".format(
                i.bases_covered, i.bases, i.percent_covered, minimum_depth, i.depth_covered))

        f.write("\nCoverage by Gene (depth >= {0})\n".format(minimum_depth))
        try:
            info = coverage.info(targeted_exons.by_transcript, minimum_depth)
        except AttributeError:
            info = coverage.calculate(targeted_exons, minimum_depth)
        for i in info:
            f.write("{0:65}{1: >3.0f}% covered (mean depth {2})\n".format(i.name.split(" ")[0], i.percent_covered, i.depth_covered))

        f.write("\nExons with incomplete coverage at depth >= {0}\n".format(minimum_depth))
        try:
            info = coverage.info(targeted_exons, minimum_depth)
        except AttributeError:
            info = coverage.calculate(targeted_exons, minimum_depth, exons=True)
        for i in info:
            if i.incompletely_covered:
                mean_depth = "(mean depth {0})".format(i.depth_covered) if (i.bases_covered>0) else ""
                f.write("{0:65}{1: >3.0f}% covered {2}\n".format(i.name.split(" ")[0], i.percent_covered, mean_depth))
    
        f.write("\nRanges not covered at  depth >= {0}\n".format(minimum_depth))
        for i in info:
            if i.incompletely_covered:
                name = i.name.split(" ")
                name = name[0] + " " + name[-1]
                f.write("{0:65}{1} (mean depth {2})\n".format(name, i.range_uncovered.locations_as_string, i.depth_uncovered))

        if "Variants" in panel or "Variants_Gene" in panel:
            header = "\nKnown mutations covered at depth >= {0} by gene\n".format(minimum_depth)
            try:
                all_variants = panel["Variants"].by_gene
            except KeyError:
                all_variants = panel["Variants_Gene"]
            targeted_variants = all_variants.overlapped_by(targeted_range)
            try:
                info = coverage.info(all_variants, minimum_depth)
            except AttributeError:
                info = coverage.calculate(all_variants, minimum_depth)
            for i in info:
                mean_depth = "(mean depth {0})".format(i.depth_covered) if (i.components_covered>0) else ""
                detectable = targeted_variants.subset(i.name).number_of_components
                f.write("{0}{1:65}{2: >3} of {6: >3} ({5: >3.0f}%) detectable mutations ({4: >3.0f}% of {3: >3} known mutations) covered {7}\n".format(
                        header, i.name, i.components_covered, i.components, i.percent_components_covered, i.components_covered*100/max(detectable,1), detectable, mean_depth))
                header = ""


        if "Variants" in panel or "Variants_Disease" in panel:
            header = "\nKnown mutations covered at depth >= {0} by disease\n".format(minimum_depth)
            try:
                all_variants = panel["Variants"].by_disease
            except KeyError:
                all_variants = panel["Variants_Disease"]

            targeted_variants = all_variants.overlapped_by(targeted_range)
            try:
                info = coverage.info(all_variants, minimum_depth)
            except AttributeError:
                info = coverage.calculate(all_variants, minimum_depth)
            for i in info:

                mean_depth = "(mean depth {0})".format(i.depth_covered) if (i.components_covered>0) else ""
                try:
                    detectable = targeted_variants.subset(i.name, exactmatch=True).number_of_components
                except AttributeError:
                    detectable = targeted_variants.subset(i.name).number_of_components
                f.write("{0}{1:65}{2: >3} of {6: >3} ({5: >3.0f}%) detectable mutations ({4: >3.0f}% of {3: >3} known mutations) covered {7}\n".format(
                        header, i.name, i.components_covered, i.components, i.percent_components_covered, i.components_covered*100/max(detectable,1), detectable, mean_depth))
                header = ""

        if "Variants" in panel or "Variants_Mutation" in panel:
            header = "\nDetectable mutations not covered at depth >= {0}\n".format(minimum_depth)
            body = []
            try:
                info = coverage.info(panel["Variants"].overlapped_by(targeted_range), minimum_depth, combine=False)
            except AttributeError:
                info = coverage.calculate(panel["Variants_Mutation"].overlapped_by(targeted_range), minimum_depth)
            for i in info:
                if i.incompletely_covered:
                    try:
                        body += ["{0}{1:10}{2:30}{3:40}({4}) (mean depth {5})\n".format(
                            header, i.gene, i.range_uncovered.locations_as_string, i.mutation, i.disease, i.depth_uncovered))]
                    except AttributeError:
                        body += ["{0}{1:10}{2:30}{3:40}({4}) (mean depth {5})\n".format(
                            header, i.name.split()[0], i.range_uncovered.locations_as_string, i.name.split()[1], i.diseases, i.depth_uncovered))]
            
            if body:
                f.write(header)
                for line in sorted(body):
                    f.write(line)








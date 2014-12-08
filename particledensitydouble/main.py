from __future__ import with_statement
import sys
import os.path
import time
from classes import *
import geometry
from fileIO import *
import version
import stringconv

#
# Functions
#


def evaluatedProfileLi(profileli):
    """ Return a list of synapses which were parsed and evaluated 
        w/o errors so far  
    """
    return [pro for pro in profileli if not pro.errflag]


def saveOutput(profileli, opt):
    """ Save a summary of results of evaluated profiles
    """

    def m(x, pixelwidth):
        return geometry.toMetricUnits(x, pixelwidth)

    def m2(x, pixelwidth):
        return geometry.toMetricUnits(x, pixelwidth ** 2)  # for area units...

    def m_inv(x):
        try:
            return 1 / m(1 / x)
        except (TypeError, ZeroDivisionError):
            return None

    def na(x):
        if x in (None, -1):
            return "N/A"
        else:
            return x


    def writeSessionSummary():
        with FileWriter("session.summary", opt) as f:
            f.writerow(["Number of evaluated profiles:", len(eval_proli)])
            if err_fli:
                f.writerow(["Number of non-evaluated profiles:", len(err_fli)])
            f.writerow(["Spatial resolution:", opt.spatial_resolution,
                        eval_proli[0].metric_unit])
            f.writerow(["Shell width:", opt.shell_width,
                        eval_proli[0].metric_unit])
            f.writerow(["Determine clusters:",
                        stringconv.yes_or_no(opt.determine_clusters)])
            if opt.determine_clusters:
                f.writerow(["Within-cluster distance:",
                            opt.within_cluster_dist,
                            eval_proli[0].metric_unit])
            if opt.monte_carlo_particle_type != 'none':
                f.writerow(["Monte Carlo simulations:", "yes"])
                f.writerow(["Monte Carlo simulations:", opt.monte_carlo_runs])
                f.writerow(["Monte Carlo simulation window:",
                            opt.monte_carlo_simulation_window])
            else:
                f.writerow(["Monte Carlo simulations:", "no"])
            if opt.determine_interpoint_dists:
                f.writerow(["Determine interpoint distances:", "yes"])
                f.writerow(["Interpoint distance mode:",
                            opt.interpoint_dist_mode.capitalize().replace('.', ' ')])
                f.writerow(["Determine shortest interpoint distances:",
                            stringconv.yes_or_no(opt.interpoint_shortest_dist)])
                f.writerow(["Determine lateral interpoint distances along "
                            "profile border:",
                            stringconv.yes_or_no(opt.interpoint_lateral_dist)])
            else:
                f.writerow(["Determine interpoint distances:", "no"])
            if clean_fli:
                f.writerow(["Input files processed cleanly:"])
                f.writerows([[fn] for fn in clean_fli])
            if nop_fli:
                f.writerow(["Input files processed but which generated no "
                            "particle distances:"])
                f.writerows([[fn] for fn in nop_fli])
            if warn_fli:
                f.writerow(["Input files processed but which generated "
                            "warnings (see log for details):"])
                f.writerows([[fn] for fn in warn_fli])
            if err_fli:
                f.writerow(["Input files not processed or not included in "
                            "summary (see log for details):"])
                f.writerows([[fn] for fn in err_fli])

    def writeProfileSummary():
        with FileWriter("profile.summary", opt) as f:
            f.writerow(["Perimeter",
                        "Area",
                        "Number of small particles (total)",
                        "Number of small particles in profile",
                        "Area density of small particles",
                        "Number of large particles (total)",
                        "Number of large particles in profile",
                        "Area density of large particles",
                        "Profile ID",
                        "Profile type",
                        "Input file",
                        "Comment"])
            f.writerows([[m(pro.path.perimeter(), pro.pixelwidth),
                          m2(pro.path.area(), pro.pixelwidth),
                          len(pro.spli),
                          len([p for p in pro.spli
                               if p.isWithinProfile(pro.path)]),
                          1e6 * (len([p for p in pro.spli
                                      if p.isWithinProfile(pro.path)])
                                 / m2(pro.path.area(), pro.pixelwidth)),
                          len(pro.lpli),
                          len([p for p in pro.lpli
                               if p.isWithinProfile(pro.path)]),
                          1e6 * (len([p for p in pro.lpli
                                      if p.isWithinProfile(pro.path)])
                                 / m2(pro.path.area(), pro.pixelwidth)),
                          pro.ID,
                          pro.profile_type,
                          os.path.basename(pro.inputfn),
                          pro.comment]
                         for pro in eval_proli])


    def writePointSummary(pType):
        if pType == "random":
            if not opt.use_random:
                return
            else:
                pli = "randomli"
                pstr = "point"
        else:
            pli = pType[0] + "pli"
            pstr = "particle"
        with FileWriter("%s.summary" % pType, opt) as f:
            f.writerow(["%s number (as appearing in input file)"
                        % pstr.capitalize(),
                        "Distance to profile border",
                        "Within profile",
                        "Profile border-associated",
                        "Profile-associated",
                        "Profile ID",
                        "Profile type",
                        "Input file",
                        "Comment"])
            f.writerows([[n + 1,
                          m(p.distToPath, pro.pixelwidth),
                          stringconv.yes_or_no(p.isWithinProfile(pro.path)),
                          stringconv.yes_or_no(p.isAssociatedWithPath),
                          stringconv.yes_or_no(p.isWithinProfile(pro.path) or
                                               p.isAssociatedWithPath),
                          pro.ID,
                          pro.profile_type,
                          os.path.basename(pro.inputfn),
                          pro.comment]
                         for pro in eval_proli
                         for n, p in enumerate(pro.__dict__[pli])])


    def writeClusterSummary(pType):
        if not opt.determine_clusters:
            return
        if pType == "small":
            cli = "spClusters"
        elif pType == "large":
            cli = "lpClusters"
        else:
            return
        with FileWriter("cluster.%s.summary" % pType, opt) as f:
            f.writerow(["Cluster number",
                        "Number of particles in cluster",
                        "Distance to profile border of centroid",
                        "Distance to nearest cluster along border",
                        "Profile ID",
                        "Profile type",
                        "Input file",
                        "Comment"])
            f.writerows([[n + 1,
                          len(c),
                          m(c.distToPath, pro.pixelwidth),
                          m(na(c.distToNearestCluster), pro.pixelwidth),
                          pro.ID,
                          pro.profile_type,
                          os.path.basename(pro.inputfn),
                          pro.comment]
                         for pro in eval_proli
                         for n, c in enumerate(pro.__dict__[cli])])


    def writeInterpointSummaries():

        def _m(x):
            return m(x, pro.pixelwidth)

        if not opt.determine_interpoint_dists:
            return
        ipRels = dict([(key, val)
                       for key, val in opt.interpoint_relations.items()
                       if val and "simulated" not in key])
        if (len(ipRels) == 0 or not
                (opt.interpoint_shortest_dist or opt.interpoint_lateral_dist)):
            return
        if not opt.use_random:
            for key, val in opt.interpoint_relations.items():
                if "random" in key and val:
                    del ipRels[key]
        if opt.interpoint_dist_mode == 'all':
            s = "all distances"
        else:
            s = "nearest neighbour distances"
        table = ["Mode: " + s]
        headerli = ipRels.keys()
        prefixli = []
        for key, val in ipRels.items():
            prefix = key[0] + key[key.index("- ") + 2] + "_"
            if prefix[0] == prefix[1]:
                prefix = prefix.replace(prefix[0], "", 1)
            prefixli.append(prefix)
        if opt.interpoint_shortest_dist and opt.interpoint_lateral_dist:
            headerli.extend(headerli)
            prefixli.extend(map(lambda s: s + "lat", prefixli))
        topheaderli = []
        if opt.interpoint_shortest_dist:
            topheaderli.append("Shortest distances")
            if opt.interpoint_lateral_dist:
                topheaderli.extend([""] * (len(ipRels) - 1))
        if opt.interpoint_lateral_dist:
            topheaderli.append("Lateral distances along profile border")
        table.extend([topheaderli, headerli])
        cols = [[] for c in prefixli]
        for pro in eval_proli:
            for n, li in enumerate([pro.__dict__[prefix + "distli"]
                                    for prefix in prefixli]):
                cols[n].extend(map(_m, li))
        # transpose cols and append to table
        table.extend(map(lambda *col: [e if e is not None else "" for e in col],
                         *cols))
        with FileWriter("interpoint.summary", opt) as f:
            f.writerows(table)


    def writeMonteCarloDistsToBorder(pType):

        def m_li(*li):
            return [m(x, pro.pixelwidth) for x in li]

        if (opt.monte_carlo_particle_type != 'both' and
                    pType not in opt.monte_carlo_particle_type):
            return
        mcli = pType[0] + "_mcli"
        table = ["Run %d" % (n + 1) for n in range(0, opt.monte_carlo_runs)]
        for pro in eval_proli:
            table.extend(map(m_li, *[[p.distToPath for p in li["pli"]]
                                     for li in pro.__dict__[mcli]]))
        with FileWriter("simulated.{}.border.distance.summary".format(pType), opt) as f:
            f.writerows(table)


    def writeMonteCarloSameIPDists(pType, distType):

        def m_li(*li):
            return [m(x, pro.pixelwidth) for x in li]

        if (opt.monte_carlo_particle_type != 'both' and
                    pType not in opt.monte_carlo_particle_type):
            return
        if (not opt.interpoint_relations["simulated %s - simulated %s"
            % (pType, pType)]):
            return
        if ((distType == "shortest" and not opt.interpoint_shortest_dist) or
                (distType == "lateral" and not opt.interpoint_lateral_dist)):
            return
        if distType == "lateral":
            shortDistType = "lat"
        else:
            shortDistType = ""
        mcli = pType[0] + "_mcli"
        table = ["Run %d" % (n + 1) for n in range(0, opt.monte_carlo_runs)]
        for pro in eval_proli:
            table.extend(map(m_li,
                             *[p for li in pro.__dict__[mcli]
                               for p in li["same"]["%sdist" % shortDistType]]))
        with FileWriter("simulated.{}.interpoint.{}.distance.summary"
                        .format(pType, distType), opt) as f:
            f.writerows(table)


    def writeMonteCarloOtherIPDists(sim_pType, distType):

        def m_li(*li):
            return [m(x, pro.pixelwidth) for x in li]

        if (opt.monte_carlo_particle_type != 'both' and
                    sim_pType not in opt.monte_carlo_particle_type):
            return
        if ((distType == "shortest" and not opt.interpoint_shortest_dist) or
                (distType == "lateral" and not opt.interpoint_lateral_dist)):
            return
        if sim_pType == "small":
            real_pType = "large"
        else:
            real_pType = "small"
        if distType == "lateral":
            shortDistType = "lat"
        else:
            shortDistType = ""
        mcli = sim_pType[0] + "_mcli"
        if opt.interpoint_relations["simulated %s - %s"
                                    % (sim_pType, real_pType)]:
            table = ["Run %d" % (n + 1) for n in range(0, opt.monte_carlo_runs)]
            for pro in eval_proli:
                table.extend(map(m_li,
                                 *[p for li in pro.__dict__[mcli]
                                   for p in li["sim-%s" % real_pType]["%sdist" % shortDistType]]))
            with FileWriter("simulated-%s.vs.%s.interpoint.%s.distance.summary"
                                    % (sim_pType, real_pType, distType), opt) as f:
                f.writerows(table)
        if opt.interpoint_relations["%s - simulated %s" % (real_pType, sim_pType)]:
            table = ["Run %d" % (n + 1) for n in range(0, opt.monte_carlo_runs)]
            for pro in eval_proli:
                table.extend(map(m_li,
                                 *[p for li in pro.__dict__[mcli]
                                   for p in li["%s-sim" % real_pType]["%sdist" % shortDistType]]))
            with FileWriter("%s.vs.simulated-%s.interpoint.%s.distance.summary"
                                    % (real_pType, sim_pType, distType), opt) as f:
                f.writerows(table)


    def writeMonteCarloClusterSummary(pType):

        if (not opt.determine_clusters or
                (opt.monte_carlo_particle_type != 'both' and
                         pType not in opt.monte_carlo_particle_type)):
            return
        mcli = pType[0] + "_mcli"
        table = ["N particles in cluster", "Run",
                 "Distance to profile border from centroid",
                 "Distance to nearest cluster",
                 "Profile ID", "Input file", "Comment"]
        for pro in eval_proli:
            for n in range(0, opt.monte_carlo_runs):
                for c in pro.__dict__[mcli][n]["clusterli"]:
                    table.extend([len(c), n + 1,
                                  m(c.distToPath, pro.pixelwidth),
                                  m(na(c.distToNearestCluster), pro.pixelwidth),
                                  pro.ID,
                                  os.path.basename(pro.inputfn),
                                  pro.comment])
        with FileWriter("simulated.%s.cluster.summary" % pType, opt) as f:
            f.writerows(table)


    sys.stdout.write("\nSaving summaries...\n")
    opt.save_result = {'any_saved': False, 'any_err': False}
    eval_proli = [pro for pro in profileli if not pro.errflag]
    clean_fli = [pro.inputfn for pro in profileli if not (pro.errflag or pro.warnflag)]
    warn_fli = [pro.inputfn for pro in profileli if pro.warnflag]
    err_fli = [pro.inputfn for pro in profileli if pro.errflag]
    nop_fli = [pro.inputfn for pro in profileli if not (pro.spli or pro.lpli)]
    writeSessionSummary()
    writeProfileSummary()
    writePointSummary("small")
    writePointSummary("large")
    writePointSummary("random")
    writeClusterSummary("small")
    writeClusterSummary("large")
    writeInterpointSummaries()
    writeMonteCarloDistsToBorder("small")
    writeMonteCarloClusterSummary("small")
    writeMonteCarloSameIPDists("small", "shortest")
    writeMonteCarloSameIPDists("small", "lateral")
    writeMonteCarloOtherIPDists("small", "shortest")
    writeMonteCarloOtherIPDists("small", "lateral")
    writeMonteCarloDistsToBorder("large")
    writeMonteCarloClusterSummary("large")
    writeMonteCarloSameIPDists("large", "shortest")
    writeMonteCarloSameIPDists("large", "lateral")
    writeMonteCarloOtherIPDists("large", "shortest")
    writeMonteCarloOtherIPDists("large", "lateral")
    if opt.save_result['any_err']:
        sys.stdout.write("Note: One or more summaries could not be saved.\n")
    if opt.save_result['any_saved']:
        sys.stdout.write("Done.\n")
    else:
        sys.stdout.write("No summaries saved.\n")

def showOptions(opt):
    sys.stdout.write("{} version: {} (Last modified {} {}, {})\n".format(
                      version.title, version.version, *version.date))
    sys.stdout.write("Output file format: %s\n" % opt.output_file_format)
    sys.stdout.write("Suffix of output files: %s\n" % opt.output_filename_ext)
    sys.stdout.write("Output directory: %s\n" % opt.output_dir)
    sys.stdout.write("Spatial resolution: %d metric units\n"
                     % opt.spatial_resolution)
    sys.stdout.write("Shell width: %d metric units\n"
                     % opt.shell_width)
    sys.stdout.write("Interpoint distances calculated: %s\n"
                     % stringconv.yes_or_no(opt.determine_interpoint_dists))
    if opt.determine_interpoint_dists:
        sys.stdout.write("Interpoint distance mode: %s\n"
                         % opt.interpoint_dist_mode.capitalize()
                         .replace('.', ' '))
        sys.stdout.write("Shortest interpoint distances: %s\n"
                         % stringconv.yes_or_no(opt.interpoint_shortest_dist))
        sys.stdout.write("Lateral interpoint distances: %s\n"
                         % stringconv.yes_or_no(opt.interpoint_lateral_dist))
    sys.stdout.write("Determine clusters: %s\n"
                     % stringconv.yes_or_no(opt.determine_clusters))
    if opt.determine_clusters:
        sys.stdout.write("Within-cluster distance: %s metric units\n"
                         % opt.within_cluster_dist)
    sys.stdout.write("Monte Carlo simulations: %s\n"
                     % stringconv.yes_or_no(
                     opt.monte_carlo_particle_type != 'none'))
    if opt.monte_carlo_particle_type != 'none':
        sys.stdout.write("Number of Monte Carlo runs: %s\n"
                         % opt.monte_carlo_runs)
        sys.stdout.write("Monte Carlo simulation window: %s\n"
                         % opt.monte_carlo_simulation_window)
        sys.stdout.write("Strict localization in simulation window: %s\n"
                         % yes_or_no(opt.monte_carlo_strict_location))
    sys.stdout.write("Clusters determined: %s\n"
                     % stringconv.yes_or_no(opt.determine_clusters))
    if opt.determine_clusters:
        sys.stdout.write("Within-cluster distance: %s metric units\n"
                         % opt.within_cluster_dist)

def getOutputFormat(opt):
    if opt.output_file_format == 'excel':
        try:
            import xls
        except ImportError:
            sys.stdout.write("Unable to write Excel files: resorting to csv "
                             "format.\n")
            opt.output_file_format = "csv"
    if opt.output_file_format == 'csv':
        opt.output_filename_ext = ".csv"
        opt.csv_format = {'dialect': 'excel', 'lineterminator': '\n',
                          'encoding': sys.getfilesystemencoding()}
        if opt.csv_delimiter == 'tab':
            opt.csv_format['delimiter'] = '\t'
    if opt.output_filename_date_suffix:
        import datetime

        opt.output_filename_suffix = "." + datetime.date.today().isoformat()
    if opt.output_filename_other_suffix != '':
        opt.output_filename_suffix += "." + opt.output_filename_other_suffix


def mainProc(parent, opt):
    """ Process profile data files
    """

    def removeDuplicateFilenames(fli):
        """ Remove duplicate filenames in input file list
        """
        for f in fli:
            if fli.count(f) > 1:
                sys.stdout.write("Duplicate input filename %s:\n   => "
                                 "removing first occurrence in list\n" % f)
                fli.remove(f)

    if not opt.input_file_list:
        sys.stdout.write("No input files.\n")
        return 0
    i, n = 0, 0
    profileli = []
    sys.stdout.write("--- Session started %s local time ---\n"
                     % time.ctime())
    removeDuplicateFilenames(opt.input_file_list)
    getOutputFormat(opt)
    showOptions(opt)
    while True:
        if i < len(opt.input_file_list):
            inputfn = opt.input_file_list[i]
            i += 1
        else:
            sys.stdout.write("\nNo more input files...\n")
            break
        parent.process_queue.put(("new_file", inputfn))
        profileli.append(ProfileData(inputfn, opt))
        profileli[-1].process(opt)
        if opt.stop_requested:
            sys.stdout.write("\n--- Session aborted by user %s local time ---\n"
                             % time.ctime())
            return 3
        if not profileli[-1].errflag:
            n += 1
            if profileli[-1].warnflag:
                sys.stdout.write("Warning(s) found while processing "
                                 "input file.\n")
                continue
        else:
            sys.stdout.write("Error(s) found while processing input file =>\n"
                             "  => No distances could be determined.\n")
            continue
    # no more input files
    errfli = [pro.inputfn for pro in profileli if pro.errflag]
    warnfli = [pro.inputfn for pro in profileli if pro.warnflag]
    if errfli:
        sys.stdout.write("\n%s input %s generated one or more "
                        "errors:\n"
                         % (stringconv.plurality("This", len(errfli)),
                            stringconv.plurality("file", len(errfli))))
        sys.stdout.write("%s\n" % "\n".join([fn for fn in errfli]))
    if warnfli:
        sys.stdout.write("\n%s input %s generated one or more warnings:\n"
                         % (stringconv.plurality("This", len(warnfli)),
                            stringconv.plurality("file", len(warnfli))))
        sys.stdout.write("%s\n" % "\n".join([fn for fn in warnfli]))
    if n > 0:
        parent.process_queue.put(("saving_summaries", ""))
        saveOutput(profileli, opt)
    else:
        sys.stdout.write("\nNo files processed.\n")
    sys.stdout.write("--- Session ended %s local time ---\n" % time.ctime())
    parent.process_queue.put(("done",""))
    opt.reset()
    if errfli:
        return 0
    elif warnfli:
        return 2
    else:
        return 1


import itertools
import os.path
import time
from .core import *
from . import geometry
from . import file_io
from . import version
from . import stringconv


#
# Functions
#

def evaluated_profile_li(profileli):
    """ Return a list of synapses which were parsed and evaluated 
        w/o errors so far  
    """
    return [pro for pro in profileli if not pro.errflag]


def save_output(profileli, opt):
    """ Save a summary of results of evaluated profiles
    """
    def m(x, pixelwidth):
        return geometry.to_metric_units(x, pixelwidth)

    def m2(x, pixelwidth):
        # For area units
        return geometry.to_metric_units(x, pixelwidth**2)

    def na(x):
        if x in (None, -1):
            return "N/A"
        else:
            return x

    def write_session_summary():
        with file_io.FileWriter("session.summary", opt) as f:
            f.writerow(["%s version:" % version.title,
                        "%s (Last modified %s %s, %s)" % ((version.version,) + version.date)])
            f.writerow(["Number of evaluated profiles:", len(eval_proli)])
            if err_fli:
                f.writerow(["Number of non-evaluated profiles:", len(err_fli)])
            f.writerow(["Metric unit:", eval_proli[0].metric_unit])
            f.writerow(["Spatial resolution:", opt.spatial_resolution,
                        eval_proli[0].metric_unit])
            f.writerow(["Shell width:", opt.shell_width, eval_proli[0].metric_unit])
            f.writerow(["Interpoint distances calculated:",
                        stringconv.yes_or_no(opt.determine_interpoint_dists)])
            if opt.determine_interpoint_dists:
                f.writerow(["Interpoint distance mode:", opt.interpoint_dist_mode])
                f.writerow(["Shortest interpoint distances:",
                            stringconv.yes_or_no(opt.interpoint_shortest_dist)])
                f.writerow(["Lateral interpoint distances:",
                            stringconv.yes_or_no(opt.interpoint_lateral_dist)])
                f.writerow(["Exclude particles outside simulation window:",
                            stringconv.yes_or_no(
                                opt.interpoint_really_exclude_particles_outside_window())])
            f.writerow(["Monte Carlo simulations:", opt.monte_carlo_particle_type.capitalize()])
            if opt.monte_carlo_particle_type != 'none':
                f.writerow(["Number of Monte Carlo runs:", opt.monte_carlo_runs])
                f.writerow(["Monte Carlo simulation window:", opt.monte_carlo_simulation_window])
                f.writerow(["Strict localization in simulation window:",
                            stringconv.yes_or_no(opt.monte_carlo_strict_location)])
            f.writerow(["Clusters determined:", stringconv.yes_or_no(opt.determine_clusters)])
            if opt.determine_clusters:
                f.writerow(["Within-cluster distance:",
                            opt.within_cluster_dist,
                            eval_proli[0].metric_unit])
            if clean_fli:
                f.writerow(["Input files processed cleanly:"])
                f.writerows([[fn] for fn in clean_fli])
            if nop_fli:
                f.writerow(["Input files processed but which generated no particle distances:"])
                f.writerows([[fn] for fn in nop_fli])
            if warn_fli:
                f.writerow(["Input files processed but which generated "
                            "warnings (see log for details):"])
                f.writerows([[fn] for fn in warn_fli])
            if err_fli:
                f.writerow(["Input files not processed or not included in "
                            "summary (see log for details):"])
                f.writerows([[fn] for fn in err_fli])

    def write_profile_summary():
        with file_io.FileWriter("profile.summary", opt) as f:
            f.writerow(["Perimeter",
                        "Area",
                        "Feret diameter",
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
            f.writerows([[m(pro.perimeter, pro.pixelwidth),
                          m2(pro.area, pro.pixelwidth),
                          m(pro.feret, pro.pixelwidth),
                          len(pro.spli),
                          len([p for p in pro.spli if p.is_within_profile]),
                          1e6 * (len([p for p in pro.spli
                                      if p.is_within_profile])
                                 / m2(pro.area, pro.pixelwidth)),
                          len(pro.lpli),
                          len([p for p in pro.lpli if p.is_within_profile]),
                          1e6 * (len([p for p in pro.lpli
                                      if p.is_within_profile])
                                 / m2(pro.area, pro.pixelwidth)),
                          pro.id,
                          pro.profile_type,
                          os.path.basename(pro.inputfn),
                          pro.comment] for pro in eval_proli])

    def write_point_summary(ptype):
        if ptype == "random":
            if not opt.use_random:
                return
            else:
                pli = "randomli"
                pstr = "point"
        else:
            pli = ptype[0] + "pli"
            pstr = "particle"
        with file_io.FileWriter("%s.summary" % ptype, opt) as f:
            f.writerow(["%s number (as appearing in input file)"
                        % pstr.capitalize(),
                        "Distance to profile border",
                        "Within profile",
                        "Profile border-associated",
                        "Profile-associated",
                        "Profile ID",
                        "Profile ptype",
                        "Input file",
                        "Comment"])
            f.writerows([[n + 1,
                          m(p.dist_to_path, pro.pixelwidth),
                          stringconv.yes_or_no(p.is_within_profile),
                          stringconv.yes_or_no(p.is_associated_with_path),
                          stringconv.yes_or_no(p.is_within_profile or
                                       p.is_associated_with_path),
                          pro.id,
                          pro.profile_type,
                          os.path.basename(pro.inputfn),
                          pro.comment] for pro in eval_proli for n, p in
                         enumerate(pro.__dict__[pli])])

    def write_cluster_summary(ptype):
        if not opt.determine_clusters:
            return
        if ptype == "small":
            cli = "s_clusterli"
        elif ptype == "large":
            cli = "l_clusterli"
        else:
            return
        with file_io.FileWriter("cluster.%s.summary" % ptype, opt) as f:
            f.writerow(["Cluster number",
                        "Number of particles in cluster",
                        "Distance to profile border of centroid",
                        "Distance to nearest cluster along border",
                        "Profile ID",
                        "Profile ptype",
                        "Input file",
                        "Comment"])
            f.writerows([[n + 1,
                          len(c),
                          m(c.dist_to_path, pro.pixelwidth),
                          m(na(c.dist_to_nearest_cluster), pro.pixelwidth),
                          pro.id,
                          pro.profile_type,
                          os.path.basename(pro.inputfn),
                          pro.comment]
                         for pro in eval_proli
                         for n, c in enumerate(pro.__dict__[cli])])

    def write_interpoint_summaries():
        if not opt.determine_interpoint_dists:
            return
        ip_rels = dict([(key, val)
                        for key, val in opt.interpoint_relations.items()
                        if val and 'simulated' not in key])
        if not opt.use_random:
            for key, val in opt.interpoint_relations.items():
                if 'random' in key and val:
                    del ip_rels[key]
        if (len(ip_rels) == 0 or not
                (opt.interpoint_shortest_dist or opt.interpoint_lateral_dist)):
            return
        table = []
        if opt.interpoint_dist_mode == 'all':
            s = "all distances"
        else:
            s = "nearest neighbour distances"
        table.append(["Mode: " + s])
        headerli = list(ip_rels.keys())
        prefixli = []
        for key, val in ip_rels.items():
            prefix = key[0] + key[key.index("- ") + 2] + "_"
            prefixli.append(prefix)
        if opt.interpoint_shortest_dist and opt.interpoint_lateral_dist:
            headerli.extend(headerli)
            prefixli.extend([t + 'lat' for t in prefixli])
        topheaderli = []
        if opt.interpoint_shortest_dist:
            topheaderli.append("Shortest distances")
            if opt.interpoint_lateral_dist:
                topheaderli.extend([""] * (len(ip_rels) - 1))
        if opt.interpoint_lateral_dist:
            topheaderli.append("Lateral distances along profile border")
        table.extend([topheaderli, headerli])
        cols = [[] for _ in prefixli]
        for pro in eval_proli:
            for n, li in enumerate([pro.__dict__[prefix + "distli"] for prefix in prefixli]):
                cols[n].extend([m(e, pro.pixelwidth) for e in li])
        # transpose cols and append to table
        table.extend(list(itertools.zip_longest(*cols, fillvalue="")))
        with file_io.FileWriter("interpoint.distances", opt) as f:
            f.writerows(table)

    def write_mc_dist_to_path(ptype):
        if opt.monte_carlo_particle_type != 'both' and ptype not in opt.monte_carlo_particle_type:
            return
        mcli = ptype[0] + "_mcli"
        table = [["Run %d" % (n + 1) for n in range(0, opt.monte_carlo_runs)]]
        for pro in eval_proli:
            table.extend(itertools.zip_longest(*[[m(p.dist_to_path, pro.pixelwidth)
                                                  for p in li['pli']]
                                                 for li in pro.__dict__[mcli]]))
        with file_io.FileWriter("simulated.{}.border.distances".format(ptype), opt) as f:
            f.writerows(table)

    def write_mc_same_ip_dists(ptype, dist_type):
        if opt.monte_carlo_particle_type != 'both' and ptype not in opt.monte_carlo_particle_type:
            return
        if not opt.interpoint_relations["simulated %s - simulated %s" % (ptype, ptype)]:
            return
        if ((dist_type == "shortest" and not opt.interpoint_shortest_dist) or
                (dist_type == "lateral" and not opt.interpoint_lateral_dist)):
            return
        if dist_type == "lateral":
            short_dist_type = "lat"
        else:
            short_dist_type = ""
        mcli = ptype[0] + "_mcli"
        table = [["Run %d" % (n + 1) for n in range(0, opt.monte_carlo_runs)]]
        for pro in eval_proli:
            table.extend(itertools.zip_longest(*[m(p, pro.pixelwidth)
                                                 for li in pro.__dict__[mcli]
                                               for p in li["same"]["%sdist" % short_dist_type]]))
        with file_io.FileWriter("simulated.{}.interpoint.{}.distances"
                                .format(ptype, dist_type), opt) as f:
            f.writerows(table)

    def write_mc_other_ip_dists(sim_ptype, dist_type):
        if (opt.monte_carlo_particle_type != 'both' and
                sim_ptype not in opt.monte_carlo_particle_type):
            return
        if ((dist_type == "shortest" and not opt.interpoint_shortest_dist) or
                (dist_type == "lateral" and not opt.interpoint_lateral_dist)):
            return
        if sim_ptype == "small":
            real_ptype = "large"
        else:
            real_ptype = "small"
        if dist_type == "lateral":
            short_dist_type = "lat"
        else:
            short_dist_type = ""
        mcli = sim_ptype[0] + "_mcli"
        if opt.interpoint_relations["%s - simulated %s" % (real_ptype, sim_ptype)]:
            table = [["Run %d" % (n + 1) for n in range(0, opt.monte_carlo_runs)]]
            for pro in eval_proli:
                table.extend(itertools.zip_longest(*[m(p, pro.pixelwidth)
                                                     for li in pro.__dict__[mcli]
                                                     for p in li["%s-sim" % real_ptype]
                                                     ["%sdist" % short_dist_type]]))
            with file_io.FileWriter("%s.vs.simulated-%s.interpoint.%s.distance.summary"
                                    % (real_ptype, sim_ptype, dist_type), opt) as f:
                f.writerows(table)

    def write_mc_cluster_summary(ptype):
        if not (opt.determine_clusters and
                (opt.monte_carlo_particle_type != 'none' and
                 ptype not in opt.monte_carlo_particle_type)):
            return
        mcli = ptype[0] + "_mcli"
        table = [["N particles in cluster", "Run",
                  "Distance to profile border from centroid",
                  "Distance to nearest cluster",
                  "Profile ID", "Input file", "Comment"]]
        for pro in eval_proli:
            for n in range(0, opt.monte_carlo_runs):
                for c in pro.__dict__[mcli][n]["clusterli"]:
                    table.append([len(c), n + 1,
                                  m(c.dist_to_path, pro.pixelwidth),
                                  m(na(c.dist_to_nearest_cluster), pro.pixelwidth),
                                  pro.id,
                                  os.path.basename(pro.inputfn),
                                  pro.comment])
        with file_io.FileWriter("simulated.%s.cluster.summary"
                                % ptype, opt) as f:
            f.writerows(table)

    def write_coords():
        if not opt.save_coords:
            return
        for pro in eval_proli:
            table = []
            table.append("# Profile coordinates including adjusted postsynaptic densities and "
                         "Monte Carlo simulated points\n")
            table.append("# %s version %s (%s %s, %s)\n" %
                         (version.title, version.version,
                          version.date[0], version.date[1], version.date[2]))
            table.append("# Generated %s\n" % time.ctime())
            table.append("INPUT_FILE %s\n" % pro.inputfn)
            table.append("IMAGE %s\n" % pro.src_img)
            table.append("PROFILE_ID %s\n" % pro.id)
            table.append("COMMENT %s\n" % pro.comment)
            table.append("PIXELWIDTH %s %s\n" %(pro.pixelwidth, pro.metric_unit))
            table.append("PROFILE_TYPE %s\n" % pro.profile_type)
            table.append("PLASMA_MEMBRANE\n")
            for p in pro.path:
                table.append("  %s, %s\n" % (p.x, p.y))
            table.append("END\n")
            for hole in pro.holeli:
                table.append("HOLE\n")
                for p in hole:
                    table.append("  %s, %s\n" % (p.x, p.y))
                table.append("END\n")
            table.append("SMALL_PARTICLES\n")
            for p in pro.spli:
                table.append("  %s, %s\n" % (p.x, p.y))
            table.append("END\n")
            table.append("LARGE_PARTICLES\n")
            for p in pro.lpli:
                table.append("  %s, %s\n" % (p.x, p.y))
            table.append("END\n")
            table.append("RANDOM_POINTS\n")
            for p in pro.randomli:
                table.append("  %s, %s\n" % (p.x, p.y))
            table.append("END\n")
            for n, cluster in enumerate(pro.s_clusterli):
                table.append("CLUSTER_CONVEX_HULL_SMALL %d\n" % (n+1))
                for p in cluster.convex_hull:
                    table.append("  %s, %s\n" % (p.x, p.y))
                table.append("END\n")
            for n, cluster in enumerate(pro.l_clusterli):
                table.append("CLUSTER_CONVEX_HULL_LARGE %d\n" % (n+1))
                for p in cluster.convex_hull:
                    table.append("  %s, %s\n" % (p.x, p.y))
                table.append("END\n")
            for n, mc in enumerate(pro.s_mcli):
                table.append("MONTE_CARLO_SMALL RUN %d\n" % (n+1))
                for p in mc['pli']:
                    table.append("  %s, %s\n" % (p.x, p.y))
                table.append("END\n")
            for n, mc in enumerate(pro.l_mcli):
                table.append("MONTE_CARLO_LARGE RUN %d\n" % (n+1))
                for p in mc['pli']:
                    table.append("  %s, %s\n" % (p.x, p.y))
                table.append("END\n")
            coords_dir = os.path.join(opt.output_dir, "coordinate_files")
            if not os.path.isdir(coords_dir):
                os.mkdir(coords_dir)
            fn = os.path.join(coords_dir,
                              os.path.basename(pro.inputfn).rstrip('.pdd') + ".coords.pdd")
            try:
                f = open(fn, "w")
                f.writelines(table)
                f.close()
            except IOError:
                sys.stderr.write("Could not write to output file %s" % fn)
            sys.stdout.write("Saved processed coordinate files to folder '%s.'" % coords_dir)

    sys.stdout.write("\nSaving summaries to %s:\n" % opt.output_dir)
    opt.save_result = {'any_saved': False, 'any_err': False}
    eval_proli = [profile for profile in profileli if not profile.errflag]
    clean_fli = [profile.inputfn for profile in profileli
                 if not (profile.errflag or profile.warnflag)]
    warn_fli = [profile.inputfn for profile in profileli if profile.warnflag]
    err_fli = [profile.inputfn for profile in profileli if profile.errflag]
    nop_fli = [profile.inputfn for profile in profileli
               if not (profile.spli or profile.lpli)]
    write_session_summary()
    write_profile_summary()
    write_point_summary("small")
    write_point_summary("large")
    write_point_summary("random")
    write_cluster_summary("small")
    write_cluster_summary("large")
    write_interpoint_summaries()
    write_mc_dist_to_path("small")
    write_mc_cluster_summary("small")
    write_mc_same_ip_dists("small", "shortest")
    write_mc_same_ip_dists("small", "lateral")
    write_mc_other_ip_dists("small", "shortest")
    write_mc_other_ip_dists("small", "lateral")
    write_mc_dist_to_path("large")
    write_mc_cluster_summary("large")
    write_mc_same_ip_dists("large", "shortest")
    write_mc_same_ip_dists("large", "lateral")
    write_mc_other_ip_dists("large", "shortest")
    write_mc_other_ip_dists("large", "lateral")
    write_coords()
    if opt.save_result['any_err']:
        sys.stdout.write("Note: One or more summaries could not be saved.\n")
    if opt.save_result['any_saved']:
        sys.stdout.write("Done.\n")
    else:
        sys.stdout.write("No summaries saved.\n")


def reset_options(opt):
    """ Deletes certain options that should always be set anew for each run
        (each time the "Start" button is pressed)
    """
    if hasattr(opt, "metric_unit"):
        delattr(opt, "metric_unit")
    if hasattr(opt, "use_grid"):
        delattr(opt, "use_grid")
    if hasattr(opt, "use_random"):
        delattr(opt, "use_random")


def show_options(opt):
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
        sys.stdout.write("Exclude particles outside simulation window: %s\n"
                         % stringconv.yes_or_no(
                              opt.interpoint_really_exclude_particles_outside_window()))
    sys.stdout.write("Determine clusters: %s\n"
                     % stringconv.yes_or_no(opt.determine_clusters))
    if opt.determine_clusters:
        sys.stdout.write("Within-cluster distance: %s metric units\n"
                         % opt.within_cluster_dist)
    sys.stdout.write("Monte Carlo simulations: %s\n"
                     % stringconv.yes_or_no(opt.monte_carlo_particle_type != 'none'))
    if opt.monte_carlo_particle_type != 'none':
        sys.stdout.write("Number of Monte Carlo runs: %s\n"
                         % opt.monte_carlo_runs)
        sys.stdout.write("Monte Carlo simulation window: %s\n"
                         % opt.monte_carlo_simulation_window)
        if opt.monte_carlo_simulation_window == "profile":
            sys.stdout.write("Strict localization in simulation window: %s\n"
                             % stringconv.yes_or_no(opt.monte_carlo_strict_location))
    sys.stdout.write("Clusters determined: %s\n"
                     % stringconv.yes_or_no(opt.determine_clusters))
    if opt.determine_clusters:
        sys.stdout.write("Within-cluster distance: %s metric units\n"
                         % opt.within_cluster_dist)


def get_output_format(opt):
    if opt.output_file_format == 'excel':
        try:
            import openpyxl
        except ImportError:
            sys.stdout.write("Unable to write Excel files: resorting to csv format.\n")
            opt.output_file_format = "csv"
    if opt.output_file_format == 'csv':
        opt.output_filename_ext = '.csv'
        opt.csv_format = {'dialect': 'excel', 'lineterminator': '\n'}
        if opt.csv_delimiter == 'tab':
            opt.csv_format['delimiter'] = '\t'
    if opt.output_filename_date_suffix:
        from datetime import date
        opt.output_filename_suffix = "." + date.today().isoformat()
    if opt.output_filename_other_suffix != '':
        opt.output_filename_suffix += "." + opt.output_filename_other_suffix


def main_proc(parent):
    """ Process profile data files
    """
    opt = parent.opt
    if not opt.input_file_list:
        file_io.sys.stdout.write("No input files.\n")
        return 0
    i, n = 0, 0
    profileli = []
    sys.stdout.write("--- Session started %s local time ---\n" % time.ctime())
    # Remove duplicate filenames
    for f in opt.input_file_list:
        if opt.input_file_list.count(f) > 1:
            sys.stdout.write("Duplicate input filename %s:\n   => "
                             "removing first occurrence in list\n" % f)
            opt.input_file_list.remove(f)
    get_output_format(opt)
    reset_options(opt)
    show_options(opt)
    while True:
        if i < len(opt.input_file_list):
            inputfn = opt.input_file_list[i]
            i += 1
        else:
            sys.stdout.write("\nNo more input files...\n")
            break
        parent.process_queue.put(("new_file", inputfn))
        profileli.append(Profile(inputfn, opt))
        profileli[-1].process()
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
        save_output(profileli, opt)
    else:
        sys.stdout.write("\nNo files processed.\n")
    sys.stdout.write("--- Session ended %s local time ---\n" % time.ctime())
    parent.process_queue.put(("done", ""))
    opt.reset()
    if errfli:
        return 0
    elif warnfli:
        return 2
    else:
        return 1
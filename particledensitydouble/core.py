import sys
import exceptions
import os.path
import random
import unicode_csv
import geometry
import file_io
import stringconv


# Convenience functions
def dot_progress(x, linelength=120, char='.'):
    """ Simple progress indicator on sys.stdout
    """
    sys.stdout.write(char)
    if (x + 1) % linelength == 0:
        sys.stdout.write('\n')

#
# Classes
#


class Point(geometry.Point):
    def __init__(self, x=None, y=None, ptype="", profile=None):
        if isinstance(x, geometry.Point):
            geometry.Point.__init__(self, x.x, x.y)
        else:
            geometry.Point.__init__(self, x, y)
        self.profile = profile
        if self.profile is not None:
            self.opt = self.profile.opt
        else:
            self.opt = None
        self.skipped = False
        self.ptype = ptype
        self.cluster = None
        self._dist_to_path = None
        self._is_within_profile = None
        self._is_within_hole = None
        self._is_associated_with_path = None
        self._is_associated_with_profile = None
        self.nearest_neighbour_dist = None
        self.nearest_neighbour_point = geometry.Point()
        self.nearest_lateral_neighbour_dist = None
        self.nearest_lateral_neighbour_point = geometry.Point()

    def determine_stuff(self):
        if self.is_within_hole:
            if self.ptype in ("sp", "lp"):
                profile_message("Discarding particle at %s: Located "
                                "within a profile hole" % self)
            self.skipped = True
            self.profile.nskipped[self.ptype] += 1
            return
        if not self.is_within_profile:
            if not self.is_within_shell:
                if self.ptype in ("sp", "lp"):
                    profile_message("Discarding particle at %s: "
                                    "Located outside the shell" % self)
                self.skipped = True
                self.profile.nskipped[self.ptype] += 1
                return

    @property
    def dist_to_path(self):
        if self._dist_to_path is None:
            self._dist_to_path = self.perpend_dist_closed_path(
                self.profile.path)
            if not self.is_within_profile:
                self._dist_to_path = -self._dist_to_path
        return self._dist_to_path

    @property
    def is_within_hole(self):
        """Determine whether self is inside a profile hole"""
        if self._is_within_hole is None:
            for h in self.profile.path.holeli:
                if self.is_within_polygon(h):
                    self._is_within_hole = True
                else:
                    self._is_within_hole = False
        return self._is_within_hole

    @property
    def is_within_profile(self):
        """Determine whether self is inside profile, excluding holes"""
        if self._is_within_profile is None:
            if (self.is_within_polygon(self.profile.path) and not
                    self.is_within_hole):
                self._is_within_profile = True
            else:
                self._is_within_profile = False
        return self._is_within_profile

    @property
    def is_within_shell(self):
        return (not self.is_within_profile and
                abs(self.dist_to_path) < geometry.to_pixel_units(
                    self.profile.opt.shell_width, self.profile.pixelwidth))

    @property
    def is_associated_with_path(self):
        if self._is_associated_with_path is None:
            if (abs(self.dist_to_path) <= geometry.to_pixel_units(
                    self.profile.opt.spatial_resolution,
                    self.profile.pixelwidth)):
                self._is_associated_with_path = True
            else:
                self._is_associated_with_path = False
        return self._is_associated_with_path

    @property
    def is_associated_with_profile(self):
        if self._is_associated_with_profile is None:
            if self.is_within_profile or self.is_associated_with_path:
                self._is_associated_with_profile = True
            else:
                self._is_associated_with_profile = False
        return self._is_associated_with_profile

    def get_nearest_neighbour(self, pointli):
        """Determine distance to nearest neighbour."""
        # Assumes that only valid (projectable, within shell etc) points
        # are in pointli
        mindist = float(sys.maxint)
        for p in pointli:
            # I will instead exclude non-desired points from the
            # supplied point list *before* calling this function
            if p is not self:
                d = self.dist(p)
                if d < mindist:
                    mindist = d
        if not mindist < float(sys.maxint):
            return None
        else:
            self.nearest_neighbour_dist = mindist
            return self.nearest_neighbour_dist

    # Todo: refactor to remove profile as argument?
    def get_nearest_lateral_neighbour(self, pointli):
        """Determine distance along profile border to nearest neighbour."""
        # Assumes that only valid (projectable, within shell etc) points
        # are in pointli
        mindist = float(sys.maxint)
        minp = Point()
        for p in pointli:
            if p is not self:
                d = self.lateral_dist_to_point(p, self.profile.path)
                if d < mindist:
                    mindist = d
                    minp = p
        if not mindist < float(sys.maxint):
            return None
        else:
            self.nearest_lateral_neighbour_dist = mindist
            self.nearest_lateral_neighbour_point = minp
            return self.nearest_lateral_neighbour_dist


class ProfileBorderData(geometry.SegmentedPath):
    def __init__(self, pointlist=None):
        if pointlist is None:
            pointlist = []
        geometry.SegmentedPath.__init__(self, pointlist)
        self.holeli = []

    def add_hole(self, pointlist=None):
        if pointlist is None:
            pointlist = []
        self.holeli.append(pointlist)
    
    def area(self):
        """Determine area of profile, excluding holes"""
        tot_hole_area = sum([h.area() for h in self.holeli])
        return geometry.SegmentedPath.area(self) - tot_hole_area

    def contains(self, p):
        """Determine if point p is inside profile, excluding holes"""
        if not p:
            return None
        return p.is_within_profile(self)


class PointList(list):
    def __init__(self, pointli, ptype, profile):
        try:
            self.extend([Point(p.x, p.y, ptype, profile) for p in pointli])
        except (AttributeError, IndexError):
            raise TypeError('not a list of Point elements')


class ClusterData(list):
    def __init__(self, pointli=None):
        if pointli is None:
            pointli = []
        try:
            self.extend([Point(p.x, p.y) for p in pointli])
        except (AttributeError, IndexError):
            raise TypeError('not a particle list')
        self.convex_hull = geometry.SegmentedPath()

    def lateral_dist_to_cluster(self, c2, border):
        """ Determine lateral distance to a cluster c2 along profile border.
        """
        path = geometry.SegmentedPath()
        c2_project, c2_seg_project = c2.convex_hull.centroid().\
            project_on_closed_path(border)
        project, seg_project = self.convex_hull.centroid().\
            project_on_closed_path(border)
        path.extend([project, c2_project])
        if c2_seg_project < seg_project:
            path.reverse()
        for n in range(min(c2_seg_project, seg_project) + 1,
                       max(c2_seg_project, seg_project)):
            path.insert(len(path)-1, border[n])
        length = path.length()
        return min(length, border.perimeter() - length)


class ProfileData:
    def __init__(self, inputfn, opt):
        self.inputfn = inputfn
        self.outputfn = ""
        self.opt = opt
        self.spli = []
        self.lpli = []
        self.randomli = []
        self.s_mcli = []
        self.l_mcli = []
        self.s_clusterli = []
        self.l_clusterli = []
        self.ss_distli, self.ss_latdistli = [], []
        self.ll_distli, self.ll_latdistli = [], []
        self.sl_distli, self.sl_latdistli = [], []
        self.ls_distli, self.ls_latdistli = [], []
        self.rs_distli, self.rs_latdistli = [], []
        self.rl_distli, self.rl_latdistli = [], []
        self.nskipped = {"sp": 0, "lp": 0, "random": 0}
        self.comment = ""
        self.pixelwidth = None
        self.metric_unit = ""
        self.posloc = geometry.Point()
        self.negloc = geometry.Point()
        self.warnflag = False
        self.errflag = False

    def process(self):
        """ Parse profile data from a file and determine distances
        """
        try:
            self._parse()
            self._check_paths()
            sys.stdout.write("Processing profile...\n")
            self.posloc = self.path.centroid()
            self.negloc = geometry.Point(None, None)
            if not self.posloc.is_within_polygon(self.path):
                self.negloc = self.posloc
                self.posloc = geometry.Point(None, None)
            for p in self.spli:                
                p.determine_stuff()
            self.spli = [p for p in self.spli if not p.skipped]
            for p in self.lpli:                
                p.determine_stuff()
            self.lpli = [p for p in self.lpli if not p.skipped]
            for r in self.randomli:
                r.determine_stuff()
            self.randomli = [r for r in self.randomli if not r.skipped]
            self._get_interdistlis()
            self._get_clusters()
            self._get_monte_carlo()
            if self.opt.stop_requested:
                return
            sys.stdout.write("Done.\n")
        except ProfileError, (self, msg):
            sys.stdout.write("Error: %s\n" % msg)
            self.errflag = True

    def _get_monte_carlo(self):
        if self.opt.monte_carlo_particle_type == "none":
            return
        sys.stdout.write("Running Monte Carlo simulations:\n")
        if self.opt.monte_carlo_particle_type in ("both", "small particles"):
            sys.stdout.write("Simulating small particles...\n")
            self.s_mcli = self._run_monte_carlo("small")
            sys.stdout.write("\n")
        if self.opt.monte_carlo_particle_type in ("both", "large particles"):
            sys.stdout.write("Simulating large particles...\n")
            self.l_mcli = self._run_monte_carlo("large")
            sys.stdout.write("\n")
            
    def _get_clusters(self):
        if not (self.opt.determine_clusters and
                self.opt.within_cluster_dist > 0):
            return
        sys.stdout.write("Determining clusters...\n")
        self.s_clusterli = self._determine_clusters(self.spli)
        self._process_clusters(self.s_clusterli)
        self.l_clusterli = self._determine_clusters(self.lpli)
        self._process_clusters(self.l_clusterli)
            
    def _get_interdistlis(self):
        if not self.opt.determine_interpoint_dists:
            return
        if True not in [val for key, val in
                        self.opt.interpoint_relations.items()
                        if "simulated" not in key]:
            return
        sys.stdout.write("Determining interpoint distances...\n")
        if self.opt.interpoint_relations["small - small"]:
            self.ss_distli, \
                self.ss_latdistli = self._get_same_interpoint_distances(
                    self.spli)

        if self.opt.interpoint_relations["large - large"]:
            self.ll_distli, \
                self.ll_latdistli = self._get_same_interpoint_distances(
                    self.lpli)
        if self.opt.interpoint_relations["small - large"]:
            self.sl_distli, \
                self.sl_latdistli = self._get_interpoint_distances2(self.spli,
                                                                    self.lpli)
        if self.opt.interpoint_relations["large - small"]:
            self.ls_distli,\
                self.ls_latdistli = self._get_interpoint_distances2(self.lpli,
                                                                    self.spli)
        if self.opt.interpoint_relations["random - small"]:
            self.rs_distli, \
                self.rs_latdistli = self._get_interpoint_distances2(
                    self.randomli, self.spli)
        if self.opt.interpoint_relations["random - large"]:
            self.rl_distli, \
                self.rl_latdistli = self._get_interpoint_distances2(
                    self.randomli, self.lpli)

    def _get_same_interpoint_distances(self, pointli):
        dli = []
        latdli = []
        for i in range(0, len(pointli)):
            if self.opt.stop_requested:
                return [], []
            if self.opt.interpoint_dist_mode == 'all':
                for j in range(i + 1, len(pointli)):
                    if self.opt.interpoint_shortest_dist:
                        dli.append(pointli[i].dist(pointli[j]))
                    if self.opt.interpoint_lateral_dist:
                        latdli.append(pointli[i].lateral_dist_to_point(
                            pointli[j], self.path))
            elif self.opt.interpoint_dist_mode == 'nearest neighbour':
                if self.opt.interpoint_shortest_dist:
                    dli.append(pointli[i].get_nearest_neighbour(pointli))
                if self.opt.interpoint_lateral_dist:
                    latdli.append(pointli[i].get_nearest_lateral_neighbour(
                        pointli))
        dli = [d for d in dli if d is not None]
        latdli = [d for d in latdli if d is not None]
        return dli, latdli

    def _get_interpoint_distances2(self, pointli, pointli2=None):
        if pointli2 is None:
            pointli2 = []
        dli = []
        latdli = []
        for i, p in enumerate(pointli):
            if self.opt.stop_requested:
                return [], []            
            if self.opt.interpoint_dist_mode == 'all':
                for p2 in pointli2:
                    if self.opt.interpoint_shortest_dist:
                        dli.append(p.dist(p2))
                    if self.opt.interpoint_lateral_dist:
                        latdli.append(p.lateral_dist_to_cluster(p2, self.path))
            elif self.opt.interpoint_dist_mode == 'nearest neighbour':
                if self.opt.interpoint_shortest_dist:
                    dli.append(p.get_nearest_neighbour(pointli2))
                if self.opt.interpoint_lateral_dist:
                    latdli.append(p.get_nearest_lateral_neighbour(pointli2))
        dli = [d for d in dli if d is not None]
        latdli = [d for d in latdli if d is not None]
        return dli, latdli

    def _run_monte_carlo(self, ptype):

        def is_valid(pt):
            if pt in mcli[n]["pli"]:
                return False
            if pt.is_within_profile:
                return True
            # The if clause below is not necessary but should speed up
            # the code a little bit
            else:
                if (self.opt.monte_carlo_simulation_window == "profile" and
                   self.opt.monte_carlo_strict_location):
                    return False
            if pt.is_within_hole:
                return False
            d = pt.perpend_dist_closed_path(self.path)
            if d is None:
                return False
            # border is set in the outer function according to
            # simulation window and opt.monte_carlo_strict_location
            return d <= border

        if self.opt.monte_carlo_simulation_window == "profile + shell":
            # Points outside shell have already been discarded
            spli = self.spli
            lpli = self.lpli
            border = geometry.to_pixel_units(self.opt.shell_width,
                                             self.pixelwidth)
        # If window == "profile"
        elif self.opt.monte_carlo_strict_location:
            spli = [p for p in self.spli if p.is_within_profile]
            lpli = [p for p in self.lpli if p.is_within_profile]
            border = 0  # just for clarity; won't actually be used
        else:
            spli = [p for p in self.spli if p.is_associated_with_profile]
            lpli = [p for p in self.lpli if p.is_associated_with_profile]
            # If shell width is smaller than spatial resolution,
            # the former must be used because all real particles
            # outside the shell have been discarded
            border = geometry.to_pixel_units(min(self.opt.shell_width,
                                                 self.opt.spatial_resolution),
                                             self.pixelwidth)
        if ptype == "small":
            numpoints = len(spli)
        elif ptype == "large":
            numpoints = len(lpli)
        else:
            return []
        box = self.path.bounding_box()
        mcli = []
        for n in range(0, self.opt.monte_carlo_runs):
            if self.opt.stop_requested:
                return []            
            dot_progress(n)            
            mcli.append({"pli": [],
                         "same": {"dist": [], "latdist": []},
                         "sim-small": {"dist": [], "latdist": []},
                         "sim-large": {"dist": [], "latdist": []},
                         "small-sim": {"dist": [], "latdist": []},
                         "large-sim": {"dist": [], "latdist": []},
                         "clusterli": []})
            p = Point()
            for i in range(0, numpoints):
                while True:
                    x = random.randint(int(box[0].x - border),
                                       int(box[1].x + border) + 1)
                    y = random.randint(int(box[0].y - border),
                                       int(box[2].y + border) + 1)
                    p = Point(x, y, "sim", self)
                    # escape the while loop when a valid simulated
                    # point is found
                    if is_valid(p):
                        break
                mcli[n]["pli"].append(p)
            for p in mcli[n]["pli"]:
                p.determine_stuff()
            if self.opt.determine_interpoint_dists:
                if (self.opt.interpoint_relations["simulated %s - simulated %s"
                                                  % (ptype, ptype)]):
                        distlis = self._get_same_interpoint_distances(
                            mcli[n]["pli"])
                        mcli[n]["same"]["dist"].append(distlis[0])
                        mcli[n]["same"]["latdist"].append(distlis[1])
                if ptype == "small":
                    ptype2 = "large"
                    pli2 = lpli
                else:
                    ptype2 = "small"
                    pli2 = spli
                if (self.opt.interpoint_relations["simulated %s - %s"
                                                  % (ptype, ptype2)]):
                    distlis = self._get_interpoint_distances2(mcli[n]["pli"],
                                                              pli2)
                    mcli[n]["sim-%s" % ptype2]["dist"].append(distlis[0])
                    mcli[n]["sim-%s" % ptype2]["latdist"].append(distlis[1])
                if (self.opt.interpoint_relations["%s - simulated %s"
                                                  % (ptype2, ptype)]):
                    distlis = self._get_interpoint_distances2(pli2,
                                                              mcli[n]["pli"])
                    mcli[n]["%s-sim" % ptype2]["dist"].append(distlis[0])
                    mcli[n]["%s-sim" % ptype2]["latdist"].append(distlis[1])
        if self.opt.determine_clusters:
            for n, li in enumerate(mcli):
                mcli[n]["clusterli"] = self._determine_clusters(li["pli"])
                self._process_clusters(mcli[n]["clusterli"])
        return mcli

    def _process_clusters(self, clusterli):
        for c in clusterli:
            if self.opt.stop_requested:
                return                
            c.convex_hull = geometry.convex_hull(c)
            hull_centroid = c.convex_hull.centroid()
            c.dist_to_path = hull_centroid.perpend_dist_closed_path(self.path)
        for c in clusterli:
            if self.opt.stop_requested:
                return            
            c.nearest_cluster = ClusterData()
            if len(clusterli) == 1:
                c.dist_to_nearest_cluster = -1
                return
            c.dist_to_nearest_cluster = sys.maxint
            for c2 in clusterli:
                if c2 != c:
                    d = c.lateral_dist_to_cluster(c2, self.path)
                    if d < c.dist_to_nearest_cluster:
                        c.dist_to_nearest_cluster = d
                        c.nearest_cluster = c2

    def _determine_clusters(self, pointli):
        """ Partition pointli into clusters; each cluster contains all points
            that are less than opt.within_cluster_dist from at least one
            other point in the cluster
        """
        clusterli = []
        for p1 in pointli:
            if self.opt.stop_requested:
                return []            
            if p1.cluster:
                continue
            for p2 in pointli:
                if p1 != p2 and p1.dist(p2) <= geometry.to_pixel_units(
                        self.opt.within_cluster_dist,
                        self.pixelwidth):
                    if p2.cluster is not None:
                        p1.cluster = p2.cluster
                        clusterli[p1.cluster].append(p1)
                        break
            else:
                p1.cluster = len(clusterli)
                clusterli.append(ClusterData([p1]))
        return clusterli

    def _parse(self):
        """Parse profile data from input file"""
        sys.stdout.write("\nParsing '%s':\n" % self.inputfn)
        li = file_io.read_file(self.inputfn)
        if not li:
            raise ProfileError(self, "Could not open input file")
        while li:
            s = li.pop(0).replace("\n", "").strip()
            if s.split(" ")[0].upper() == "IMAGE":
                self.src_img = s.split(" ")[1]
            elif s.split(" ")[0].upper() == "PROFILE_ID":
                try:
                    self.ID = s.split(" ")[1]
                except IndexError:
                    self.ID = ''
            elif s.split(" ")[0].upper() == "COMMENT":
                try:
                    self.comment = s.split(" ", 1)[1]
                except IndexError:
                    self.comment = ''
            elif s.split(" ")[0].upper() == "PROFILE_TYPE":
                try:
                    self.profile_type = s.split(" ", 1)[1]
                except IndexError:
                    self.profile_type = ''                    
            elif s.split(" ")[0].upper() == "PIXELWIDTH":
                try: 
                    self.pixelwidth = float(s.split(" ")[1])
                    self.metric_unit = s.split(" ")[2]
                except (IndexError, ValueError):
                    raise ProfileError(self, 
                                       "PIXELWIDTH is not a valid number")
            elif s.upper() == "PROFILE_BORDER":
                self.path = ProfileBorderData(self._get_coords(li, "path"))
            elif s.upper() in ("PROFILE_HOLE", "HOLE"):
                self.path.add_hole(geometry.SegmentedPath(
                    self._get_coords(li, "hole")))
            elif s.upper() == "SMALL_PARTICLES":
                self.spli = PointList(self._get_coords(li, "sp"), "sp", self)
            elif s.upper() == "LARGE_PARTICLES":
                self.lpli = PointList(self._get_coords(li, "lp"), "lp", self)
            elif s.upper() in ("RANDOM_POINTS", "RANDOM"):
                self.randomli = PointList(self._get_coords(li, "random"),
                                          "random", self)
            elif s[0] != "#":          # unless specifically commented out
                profile_warning(self, "Unrecognized string '" + s +
                                      "' in input file")
        # Now, let's see if everything was found
        self._check_parsed_data()

    def _check_parsed_data(self):
        """See if the synapse data was parsed correctly, and print info
        on the parsed data to stdout.
        """
        self._check_var_default(self, 'src_img', "Source image", "N/A")
        self._check_var_default(self, 'ID', "Profile ID", "N/A")
        self._check_var_default(self, 'comment', "Comment", "")
        self._check_var_val(self, 'metric_unit', "Metric unit", 'metric_unit')
        self._check_required_var(self, 'pixelwidth', "Pixel width",
                                 self.metric_unit)
        self._check_var_default(self, 'profile_type', "Profile type", "N/A")
        self._check_list_var(self, 'path', 'Profile border', 'nodes', 2)
        self._check_list_var(self, 'spli', 'Small points', '', 0)
        self._check_list_var(self, 'lpli', 'Large points', '', 0)
        self._check_table_var(self.path, 'holeli', "Hole", "Holes", 0, 2)
        self._check_var_exists(self, 'randomli', "Random points", 'use_random')

    def _check_required_var(self, parent, var_to_check, var_str, post_str):
        """Confirm that parent has a required variable; else, raise
        ProfileError.
        """
        if not hasattr(parent, var_to_check):
            raise ProfileError(self, "%s not found in input file" % var_str)
        else:
            sys.stdout.write("  %s: %s %s\n"
                             % (var_str, parent.__dict__[var_to_check],
                                post_str))

    @staticmethod
    def _check_list_len(var, min_len):
        """Return True if var is a list and has at least min_len
        elements, else False.
        """
        return isinstance(var, list) and len(var) >= min_len

    def _check_list_var(self, parent, var_to_check, var_str, post_str, min_len):
        """Confirms that parent has a var_to_check that is a list and
        has at least min_len elements; if var_to_check does not exist
        and min_len <= 0, assigns an empty list to var_to_check. Else,
        raise a ProfileError.
        """
        if not hasattr(parent, var_to_check):
            if min_len > 0:
                raise ProfileError(self, "%s not found in input file"
                                   % var_str)
            else:
                parent.__dict__[var_to_check] = []
        elif not self._check_list_len(parent.__dict__[var_to_check], min_len):
            raise ProfileError(self, "%s has too few coordinates" % var_str)
        if post_str != '':
            post_str = ' ' + post_str
        sys.stdout.write("  %s%s: %d\n"
                         % (var_str, post_str,
                            len(parent.__dict__[var_to_check])))

    def _check_table_var(self, parent, var_to_check, var_str_singular,
                         var_str_plural, min_len_1, min_len_2):
        """Confirms that var_to_check exists, is a list and has at
        least min_len_1 elements, and that each of these has at least
        min_len_2 subelements; if var_to_check does not exist and
        min_len_1 <= 0, assigns an empty list to var_to_check. Else,
        raise ProfileError.
        """
        if not hasattr(parent, var_to_check):
            if min_len_1 > 0:
                raise ProfileError(self, "%s not found in input file"
                                   % var_str_plural)
            else:
                parent.__dict__[var_to_check] = []
        elif not self._check_list_len(parent.__dict__[var_to_check], min_len_1):
            raise ProfileError(self, "Too few %s found in input file"
                               % var_str_plural.lower())
        else:
            for element in parent.__dict__[var_to_check]:
                if not self._check_list_len(element, min_len_2):
                    raise ProfileError(self, "%s has too few coordinates"
                                       % var_str_singular)
        sys.stdout.write("  %s: %d\n" % (var_str_plural,
                                         len(parent.__dict__[var_to_check])))

    @staticmethod
    def _check_var_default(parent, var_to_check, var_str, default=""):
        """Checks if var_to_check exists; if not, assign the default
        value to var_to_check. Never raises a ProfileError.
        """
        if not hasattr(parent, var_to_check):
            parent.__dict__[var_to_check] = default
        sys.stdout.write("  %s: %s\n" % (var_str,
                                         parent.__dict__[var_to_check]))

    def _check_var_exists(self, parent, var_to_check, var_str, optflag):
        """Checks for consistency between profiles with respect to the
        existence of var_to_check (i.e., var_to_check must be present
        either in all profiles or in none).

        If optflag is not set (i.e., this is the first profile), then
        set optflag to True or False depending on the existence of
        var_to_check. If optflag is already set (for consequent
        profiles), var_to_check must (if optflag is True) or must not
        (if optflag is False) exist. If not so, raise ProfileError.
        """
        if not hasattr(parent.opt, optflag):
            if hasattr(self, var_to_check):
                parent.opt.__dict__[optflag] = True
            else:
                parent.opt.__dict__[optflag] = False
        if parent.opt.__dict__[optflag]:
            if hasattr(parent, var_to_check):
                sys.stdout.write("  %s: yes\n" % var_str)
            else:
                raise ProfileError(self, "%s not found in input file" % var_str)
        elif hasattr(parent, var_to_check):
            raise ProfileError(self, "%s found but not expected" % var_str)
        else:
            sys.stdout.write("  %s: no\n" % var_str)

    def _check_var_val(self, parent, var_to_check, var_str, optvar):
        """Checks for consistency between profiles with respect to the
        value of var_to_check (i.e., var_to_check must be present and
        have equal value in all profiles).

        If optvar is not set (i.e., this is the first profile), then
        set optflag to the value of var_to_check. If optvar is already
        set (for consequent profiles), the value of var_to_check must
        be equal to that of optvar. If not so, raise ProfileError.
        """
        if not hasattr(parent, var_to_check):
            raise ProfileError(self, "%s not found in input file" % var_str)
        if not hasattr(parent.opt, optvar):
            parent.opt.__dict__[optvar] = parent.__dict__[var_to_check]
        elif parent.__dict__[var_to_check] == parent.opt.__dict__[optvar]:
            pass  # really no point in pointing out that it's ok
            # sys.stdout.write("  %s: %s\n"
            #                  % (var_str, parent.__dict__[var_to_check]))
        else:
            raise ProfileError(self, "%s value '%s'  differs from the value "
                                     "specified ('%s') in the first input file"
                               % (var_str, parent.__dict__[var_to_check],
                                  parent.opt.__dict__[optvar]))

    def _check_paths(self):
        """Check if profile border and holes intersect with themselves."""

        def _check_path(_path, s):
            for p in range(0, len(_path) - 3):
                for q in range(0, len(_path) - 1):
                    if p not in (q, q + 1) and p + 1 not in (q, q + 1):
                        if geometry.segment_intersection(_path[p],
                                                         _path[p + 1],
                                                         _path[q],
                                                         _path[q + 1]):
                            raise ProfileError(
                                self, "%s invalid (crosses itself)" % s)
            return True

        _check_path(self.path, "Profile border")
        for path in self.path.holeli:
            _check_path(path, "Hole")
        for n, h in enumerate(self.path.holeli):
            if not h.is_simple_polygon():
                raise ProfileError(self,
                                   "Profile hole %d is not a simple polygon"
                                   % (n + 1))
            for n2, h2 in enumerate(self.path.holeli[n + 1:]):
                if h.overlaps_polygon(h2):
                    raise ProfileError(self,
                                       "Profile hole %d overlaps with hole %d "
                                       % (n + 1, n + n2 + 2))
        sys.stdout.write("  Paths are ok.\n")

    def _get_coords(self, strli, coord_type=""):
        """ Pop point coordinates from list strli.
            When an element of strli is not a valid point,
            a warning is issued.
        """
        pointli = []
        s = strli.pop(0).replace("\n", "").replace(" ", "").strip()
        while s != "END":
            try:
                p = geometry.Point(float(s.split(",")[0]),
                                   float(s.split(",")[1]))
                if pointli and (p == pointli[-1] or
                                (coord_type in ('sp', 'lp', 'random')
                                 and p in pointli)):
                    sys.stdout.write("Duplicate %s coordinates %s: skipping "
                                     "2nd instance\n" % (coord_type, p))
                else:
                    pointli.append(p)                    
            except ValueError:
                if s[0] != "#":
                    profile_warning(self, "'%s' not valid %s coordinates"
                                    % (s, coord_type))
                else:
                    pass 
            s = strli.pop(0).replace("\n", "").strip()
        # For some reason, sometimes the endnodes have the same coordinates;
        # in that case, delete the last endnode to avoid division by zero
        if (len(pointli) > 1) and (pointli[0] == pointli[-1]):
            del pointli[-1]
        return pointli        

# end of class Profile
        

class OptionData:
    def __init__(self):
        self.input_file_list = []
        self.spatial_resolution = 25
        self.shell_width = 0   # Skip points farther than this from profile
        self.outputs = {'profile summary': True, 'point summary': True,
                        'random summary': True, 'session summary': True}
        self.output_file_format = 'excel'
        self.output_filename_ext = '.xls'
        self.input_filename_ext = ('.pd', '.pdd')
        self.output_filename_suffix = ''
        self.output_filename_other_suffix = ''
        self.output_filename_date_suffix = True
        self.output_filename_use_other_suffix = False
        self.csv_delimiter = 'comma'
        self.action_if_output_file_exists = 'overwrite'
        self.output_dir = ''
        self.determine_clusters = False
        self.within_cluster_dist = 100
        self.monte_carlo_particle_type = 'none'
        self.monte_carlo_runs = 99
        self.monte_carlo_simulation_window = 'profile'
        self.monte_carlo_strict_location = False
        self.determine_interpoint_dists = False
        self.interpoint_dist_mode = 'nearest neighbour'
        self.interpoint_relations = {'small - small': True,
                                     'large - large': True,
                                     'small - large': True,
                                     'large - small': True,
                                     'random - small': True,
                                     'random - large': True,
                                     'simulated small - simulated small': False,
                                     'simulated large - simulated large': False,
                                     'small - simulated large': False,
                                     'large - simulated small': False,
                                     'simulated small - large': False,
                                     'simulated large - small': False}
        self.interpoint_shortest_dist = True
        self.interpoint_lateral_dist = False
        self.stop_requested = False
        
    def reset(self):
        """ Resets all options to default, and removes those that are not 
            set in __init__(). 
        """
        self.__dict__ = {}
        self.__init__()
# end of class OptionData


class ProfileError(exceptions.Exception):
    def __init__(self, profile, msg):
        self.args = (profile, msg + ".")


def profile_warning(profile, msg):
    """ Issue a warning
    """
    sys.stdout.write("Warning: %s.\n" % msg)
    profile.warnflag = True      


def profile_message(msg):
    """ Show a message
    """
    sys.stdout.write("%s.\n" % msg)

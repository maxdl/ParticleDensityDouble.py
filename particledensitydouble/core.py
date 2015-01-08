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
    def __init__(self, x=None, y=None, ptype=""):
        if isinstance(x, geometry.Point):
            geometry.Point.__init__(self, x.x, x.y)
        else:
            geometry.Point.__init__(self, x, y)
        self.skipped = False
        self.ptype = ptype
        self.cluster = None
        self.dist_to_path = None
        self.__is_within_profile = None
        self.is_associated_with_path = None
        self.is_associated_with_profile = None
        self.nearest_neighbour_dist = None
        self.nearest_neighbour_point = geometry.Point()
        self.nearest_lateral_neighbour_dist = None
        self.nearest_lateral_neighbour_point = geometry.Point()

    def determine_stuff(self, profile):
        if self.__is_within_hole(profile.path):
            if self.ptype == "gold":
                profile_message("Discarding particle at %s: Located "
                                "within a profile hole" % self)
            self.skipped = True
            profile.nskipped[self.ptype] += 1
            return
        self.dist_to_path = self.perpend_dist_closed_path(profile.path)
        if not self.is_within_polygon(profile.path):
            if self.dist_to_path > geometry.to_pixel_units(
                    profile.opt.shell_width,
                    profile.pixelwidth):
                if self.ptype == "gold":
                    profile_message("Discarding point at %s: "
                                    "Located outside the shell" % self)
                self.skipped = True
                profile.nskipped[self.ptype] += 1
                return                 
            self.dist_to_path = -self.dist_to_path
        self.__is_within_profile = self.__is_within_profile(profile.path)
        self.is_associated_with_path = (abs(self.dist_to_path)
                                        <= geometry.to_pixel_units(
                                        profile.opt.spatial_resolution,
                                        profile.pixelwidth))
        self.is_associated_with_profile = (self.is_within_profile
                                           or self.is_associated_with_path)
                                                     
    def __is_within_hole(self, path):
        """Determine whether self is inside a profile hole"""
        for h in path.holeli:
            if self.is_within_polygon(h):
                return True
        return False

    def is_within_profile(self, path):
        """Determine whether self is inside profile, excluding holes"""
        if (not self.is_within_polygon(path)) or self.__is_within_hole(path):
            return False
        return True

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

    def get_nearest_lateral_neighbour(self, pointli, profile):
        """Determine distance along profile border to nearest neighbour."""
        # Assumes that only valid (projectable, within shell etc) points
        # are in pointli
        mindist = float(sys.maxint)
        minp = Point()
        for p in pointli:
            if p is not self:
                d = self.lateral_dist_to_point(p, profile.path)
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
    def __init__(self, pointli, ptype):
        try:
            self.extend([Point(p.x, p.y, ptype) for p in pointli])
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
        c2_project, c2_seg_project = c2.convexHull.centroid().\
            project_on_closed_path(border)
        project, seg_project = self.convex_hull.centroid().\
            project_on_closedPath(border)
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
            self.__parse()
            self.__check_paths()
            sys.stdout.write("Processing profile...\n")
            self.posloc = self.path.centroid()
            self.negloc = geometry.Point(None, None)
            if not self.posloc.is_within_polygon(self.path):
                self.negloc = self.posloc
                self.posloc = geometry.Point(None, None)
            for p in self.spli:                
                p.determine_stuff(self.spli, self)
            self.spli = [p for p in self.spli if not p.skipped]
            for p in self.lpli:                
                p.determine_stuff(self.lpli, self)
            self.lpli = [p for p in self.lpli if not p.skipped]
            for r in self.randomli:
                r.determine_stuff(self.randomli, self)
            self.randomli = [r for r in self.randomli if not r.skipped]
            self.get_interdistlis()
            self.get_clusters()
            self.get_monte_carlo()
            if self.opt.stop_requested:
                return
            sys.stdout.write("Done.\n")
        except ProfileError, (self, msg):
            sys.stdout.write("Error: %s\n" % msg)
            self.errflag = True

    def get_monte_carlo(self):
        if self.opt.monte_carlo_particle_type == "none":
            return
        sys.stdout.write("Running Monte Carlo simulations:\n")
        if self.opt.monte_carlo_particle_type in ("both", "small particles"):
            sys.stdout.write("Simulating small particles...\n")
            self.s_mcli = self.__run_monte_carlo("small")
            sys.stdout.write("\n")
        if self.opt.monte_carlo_particle_type in ("both", "large particles"):
            sys.stdout.write("Simulating large particles...\n")
            self.l_mcli = self.__run_monte_carlo("large")
            sys.stdout.write("\n")
            
    def get_clusters(self):
        if not (self.opt.determine_clusters and
                self.opt.within_cluster_dist > 0):
            return
        sys.stdout.write("Determining clusters...\n")
        self.s_clusterli = self.determine_clusters(self.spli)
        self.process_clusters(self.s_clusterli)
        self.l_clusterli = self.determine_clusters(self.lpli)
        self.process_clusters(self.l_clusterli)
            
    def get_interdistlis(self):
        if not self.opt.determine_interpoint_dists:
            return
        if True not in [val for key, val in
                        self.opt.interpoint_relations.items()
                        if "simulated" not in key]:
            return
        sys.stdout.write("Determining interpoint distances...\n")
        if self.opt.interpoint_relations["small - small"]:
            self.ss_distli, \
                self.ss_latdistli = self.get_same_interpoint_distances(
                    self.spli)

        if self.opt.interpoint_relations["large - large"]:
            self.ll_distli, \
                self.ll_latdistli = self.get_same_interpoint_distances(
                    self.lpli)
        if self.opt.interpoint_relations["small - large"]:
            self.sl_distli, \
                self.sl_latdistli = self.get_interpoint_distances2(self.spli,
                                                                   self.lpli)
        if self.opt.interpoint_relations["large - small"]:
            self.ls_distli,\
                self.ls_latdistli = self.get_interpoint_distances2(self.lpli,
                                                                   self.spli)
        if self.opt.interpoint_relations["random - small"]:
            self.rs_distli, \
                self.rs_latdistli = self.get_interpoint_distances2(
                    self.randomli, self.spli)
        if self.opt.interpoint_relations["random - large"]:
            self.rl_distli, \
                self.rl_latdistli = self.get_interpoint_distances2(
                    self.randomli, self.lpli)

    def get_same_interpoint_distances(self, pointli):
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
                        latdli.append(pointli[i].lateral_dist_to_cluster(
                            pointli[j], self.path))
            elif self.opt.interpoint_dist_mode == 'nearest neighbour':
                if self.opt.interpoint_shortest_dist:
                    dli.append(pointli[i].get_nearest_neighbour(pointli, self))
                if self.opt.interpoint_lateral_dist:
                    latdli.append(pointli[i].get_nearest_lateral_neighbour(
                        pointli, self))
        dli = [d for d in dli if d is not None]
        latdli = [d for d in latdli if d is not None]
        return dli, latdli

    def get_interpoint_distances2(self, pointli, pointli2=None):
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
                    dli.append(p.get_nearest_neighbour(pointli2, self))
                if self.opt.interpoint_lateral_dist:
                    latdli.append(p.get_nearest_lateral_neighbour(pointli2,
                                                                  self))
        dli = [d for d in dli if d is not None]
        latdli = [d for d in latdli if d is not None]
        return dli, latdli

    def __run_monte_carlo(self, ptype):

        def is_valid(pt):
            if pt in mcli[n]["pli"]:
                return False
            if pt.is_within_profile(self.path):
                return True
            if pt.__is_within_hole(self.path):
                return False
            d = pt.perpend_dist_closed_path(self.path)
            if d is None:
                return False
            if (self.opt.monte_carlo_simulation_window == "profile" and not
                self.opt.monte_carlo_strict_location and
                (abs(d) < geometry.to_pixel_units(self.opt.spatial_resolution,
                                                  self.pixelwidth))):
                    return True
            if (self.opt.monte_carlo_simulation_window == "profile + shell" and
                    abs(d) <= geometry.to_pixel_units(self.opt.shell_width,
                                                      self.pixelwidth)):
                return True
            return False

        # TODO: check if correct if using monte_carlo_strict_location or not
        if self.opt.monte_carlo_simulation_window == "profile + shell":
            spli = self.spli    # Points outside the shell have already
            lpli = self.lpli    # been skipped.
        else:                   # If window == "profile"
            spli = [p for p in self.spli if p.is_within_profile]
            lpli = [p for p in self.lpli if p.is_within_profile]
        if ptype == "small":
            numpoints = len(spli)
        elif ptype == "large":
            numpoints = len(lpli)
        else:
            return []
        box = self.path.bounding_box()
        border = 0
        if self.opt.monte_carlo_simulation_window == "profile + shell":
            border = geometry.to_pixel_units(self.opt.shell_width,
                                             self.pixelwidth)
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
                    p = Point(x, y)
                    # escape the while loop when a valid simulated
                    # point is found
                    if is_valid(p):
                        break
                mcli[n]["pli"].append(p)
            for p in mcli[n]["pli"]:
                p.determine_stuff(mcli[n]["pli"], self)
            if self.opt.determine_interpoint_dists:
                if (self.opt.interpoint_relations["simulated %s - simulated %s"
                                                  % (ptype, ptype)]):
                        distlis = self.get_same_interpoint_distances(
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
                    distlis = self.get_interpoint_distances2(mcli[n]["pli"],
                                                             pli2)
                    mcli[n]["sim-%s" % ptype2]["dist"].append(distlis[0])
                    mcli[n]["sim-%s" % ptype2]["latdist"].append(distlis[1])
                if (self.opt.interpoint_relations["%s - simulated %s"
                                                  % (ptype2, ptype)]):
                    distlis = self.get_interpoint_distances2(pli2,
                                                             mcli[n]["pli"])
                    mcli[n]["%s-sim" % ptype2]["dist"].append(distlis[0])
                    mcli[n]["%s-sim" % ptype2]["latdist"].append(distlis[1])
        if self.opt.determine_clusters:
            for n, li in enumerate(mcli):
                mcli[n]["clusterli"] = self.determine_clusters(li["pli"])
                self.process_clusters(mcli[n]["clusterli"])
        return mcli

    def process_clusters(self, clusterli):
        for c in clusterli:
            if self.opt.stop_requested:
                return                
            c.convex_hull = geometry.convex_hull(c)
            hull_centroid = c.convex_hull.centroid()
            c.distToPath = hull_centroid.centroid().perpend_dist_closed_path(
                self.path)
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

    def determine_clusters(self, pointli):
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

    def __parse(self):
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
                self.path = ProfileBorderData(self.__get_coords(li, "path"))
            elif s.upper() == "PROFILE_HOLE":
                self.path.add_hole(geometry.SegmentedPath(
                    self.__get_coords(li, "hole")))
            elif s.upper() == "SMALL_PARTICLES":
                self.spli = PointList(self.__get_coords(li, "small particle"),
                                      "gold")
            elif s.upper() == "LARGE_PARTICLES":
                self.lpli = PointList(self.__get_coords(li, "large particle"),
                                      "gold")
            elif s.upper() == "RANDOM_POINTS":
                self.randomli = PointList(self.__get_coords(li, "random"),
                                          "random")
            elif s[0] != "#":          # unless specifically commented out
                profile_warning(self, "Unrecognized string '" + s +
                                      "' in input file")
        # Now, let's see if everything was found
        self.check_parsed_data()

    def check_parsed_data(self):
        """See if the synapse data was parsed correctly, and print info
        on the parsed data to stdout.
        """
        self.check_var_default(self, 'src_img', "Source image", "N/A")
        self.check_var_default(self, 'ID', "Profile ID", "N/A")
        self.check_var_default(self, 'comment', "Comment", "")
        self.check_var_val(self, 'metric_unit', "Metric unit", 'metric_unit')
        self.check_required_var(self, 'pixelwidth', "Pixel width",
                                self.metric_unit)
        self.check_var_default(self, 'profile_type', "Profile type", "N/A")
        self.check_list_var(self, 'path', 'Profile border', 'nodes', 2)
        self.check_list_var(self, 'spli', 'Small points', '', 0)
        self.check_list_var(self, 'lpli', 'Large points', '', 0)
        self.check_table_var(self.path, 'holeli', "Hole", "Holes", 0, 2)
        self.check_var_exists(self, 'randomli', "Random points", 'use_random')

    def check_required_var(self, parent, var_to_check, var_str, post_str):
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
    def check_list_len(var, min_len):
        """Return True if var is a list and has at least min_len
        elements, else False.
        """
        return isinstance(var, list) and len(var) >= min_len

    def check_list_var(self, parent, var_to_check, var_str, post_str, min_len):
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
        elif not self.check_list_len(parent.__dict__[var_to_check], min_len):
            raise ProfileError(self, "%s has too few coordinates" % var_str)
        if post_str != '':
            post_str = " " + post_str
        sys.stdout.write("  %s%s: %d\n"
                         % (var_str, post_str,
                            len(parent.__dict__[var_to_check])))

    def check_table_var(self, parent, var_to_check, var_str_singular,
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
        elif not self.check_list_len(parent.__dict__[var_to_check], min_len_1):
            raise ProfileError(self, "Too few %s found in input file"
                               % var_str_plural.lower())
        else:
            for element in parent.__dict__[var_to_check]:
                if not self.check_list_len(element, min_len_2):
                    raise ProfileError(self, "%s has too few coordinates"
                                       % var_str_singular)
        sys.stdout.write("  %s: %d\n" % (var_str_plural,
                                         len(parent.__dict__[var_to_check])))

    @staticmethod
    def check_var_default(parent, var_to_check, var_str, default=""):
        """Checks if var_to_check exists; if not, assign the default
        value to var_to_check. Never raises a ProfileError.
        """
        if not hasattr(parent, var_to_check):
            parent.__dict__[var_to_check] = default
        sys.stdout.write("  %s: %s\n" % (var_str,
                                         parent.__dict__[var_to_check]))

    def check_var_exists(self, parent, var_to_check, var_str, optflag):
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

    def check_var_val(self, parent, var_to_check, var_str, optvar):
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

    def __check_paths(self):
        """Check if profile border and holes intersect with themselves."""

        def check_path(_path, s):
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

        check_path(self.path, "Profile border")
        for path in self.path.holeli:
            check_path(path, "Hole")
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

    def __get_coords(self, strli, coord_type=""):
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
                                (coord_type in ('small particle',
                                                'large particle', 'random')
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
        self.input_filename_ext = '.pdd'
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

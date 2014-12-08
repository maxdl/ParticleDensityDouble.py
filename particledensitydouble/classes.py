import sys
import random
import exceptions
import os.path
import csv
from geometry import *
from fileIO import *
from stringconv import *

#
# Global constants
# 

DefaultOutputFileFormat = 'csv'
DefaultOutputFilenameExt = '.csv'
DefaultCSVDelimiter = 'comma'
DefaultSpatialResolution = 25  # metric length units 

# Convenience functions
def dotProgress(x, linelength=120, char='.'):
    """ Simple progress indicator on sys.stdout
    """
    sys.stdout.write(char)
    if (x + 1) % linelength == 0:
        sys.stdout.write('\n')

#
# Classes
#
class Particle(Point):    
    def __init__(self, x=None, y=None, type=""):
        Point.__init__(self, x, y)
        self.__skipped = False
        self.type = type
        self.cluster = None

    def skip(self, skipped):
        self.__skipped = skipped
    
    def is_skipped(self):
        return self.__skipped
    
    skipped = property(is_skipped, skip)

    def determineStuff(self, parentli, profile, opt):
        if self.isWithinHole(profile.path):
            if self.type == "gold":
                ProfileMessage(profile, "skipping point at %s: Located "
                                        "within a profile hole" % self)
            self.skipped = True
            return
        self.distToPath = self.perpendDistClosedPath(profile.path)                                        
        if not self.isWithinPolygon(profile.path):
            if self.distToPath > toPixelUnits(opt.shell_width,
                                              profile.pixelwidth):
                if self.type == "gold":
                    ProfileMessage(profile, "Skipping point at %s: "
                                            "Located farther "
                                            "than %d %s from profile" 
                                            % (self, opt.shell_width,
                                               profile.metric_unit))
                self.skipped = True
                return                 
            self.distToPath = -self.distToPath
        self.isAssociatedWithPath = (abs(self.distToPath)
                                     <= toPixelUnits(opt.spatial_resolution, 
                                                     profile.pixelwidth))
        self.isAssociatedWithProfile = (self.isWithinProfile(profile.path) or
                                        self.isAssociatedWithPath)
                                                     


    def isWithinHole(self, path):
        """  Determine whether self is inside a profile hole
        """                
        for h in path.holeli:
            if self.isWithinPolygon(h):
                return True
        return False
   
   
    def isWithinProfile(self, path):
        """  Determine whether self is inside profile, excluding holes
        """        
        if (not self.isWithinPolygon(path)) or self.isWithinHole(path):
            return False
        return True



    def determineNearestNeighbour(self, pointli, profile, opt):
        nearestNeighbourDist = None
        nearestNeighbourParticle = Particle()
        #if not self.isWithinProfile:       I will check this *before* calling
        #    return None                    this function I think
        mindist = float(sys.maxint)
        for p in pointli:
            # I will instead exclude non-desired points from the supplied point
            # list *before* calling this function
            #if p is not self and p.isWithinProfile:
            if p is not self:
                d = self.dist(p)
                if d < mindist:
                    mindist = d
                    minp = p
        if not mindist < float(sys.maxint):
            return None
        else:
            nearestNeighbourDist = mindist
            nearestNeighbourParticle = minp
            return nearestNeighbourDist

    def determineNearestLateralNeighbour(self, pointli, profile, opt):
        nearestLateralNeighbourDist = None
        nearestLateralNeighbourParticle = Particle()
        #if not self.isWithinProfile:      I will check this *before* calling
        #    return None                   this function I think
        mindist = float(sys.maxint)
        for p in pointli:
            # I will instead exclude non-desired points from the supplied point
            # list *before* calling this function
            #if p is not self and p.isWithinProfile:
            if p is not self:
                d = self.lateralDist(p, profile.path)
                if d < mindist:
                    mindist = d
                    minp = p
        if not mindist < float(sys.maxint):
            return None
        else:
            nearestLateralNeighbourDist = mindist
            nearestLateralNeighbourParticle = minp
            return nearestLateralNeighbourDist

    def perpendDist(self, m, negloc=Point(None, None), 
                    posloc=Point(None, None), 
                    doNotCareIfOnOrOffSeg=False):        
        """" Calculate distance from the point to an open path m;
             negloc is a point defined to have a negative distance
             to the path; posloc is a point defined to have a positive 
             distance to the path; if neither negloc nor posloc is 
             defined, absolute distance is returned.             
        """
        mindist = float(sys.maxint)
        on_M = False
        for n in range(0, len(m) - 1):
            if (m[n].x != -1) and (m[n+1].x != -1):
                on_this_seg, d = self.distToSegment(m, n)
                if d <= mindist:      # smallest distance so far...
                    mindist = d
                    if on_this_seg or doNotCareIfOnOrOffSeg: 
                        on_M = True   # least distance and "on" segment (not
                                      # completely true; see distToSegment())
                    else:
                        on_M = False      # least distance but "off" segment
        if not on_M:
            return None     
        # If polarity (posloc or negloc) is defined, we say that points
        # on the positive side of the path have positive distances to the 
        # path, while other points have negative distances. To
        # determine this, we count the number of path segments 
        # dissected by the line between the point and negloc (posloc).
        # Even (odd) number => the point and negloc (posloc) are on the same
        # same side of the path; odd number => different side.
        if negloc and self.segmentCrossingNumber(m, negloc) % 2 == 0:
            mindist = -mindist
        elif posloc and self.segmentCrossingNumber(m, posloc) % 2 != 0:
            mindist = -mindist
        return mindist

    def lateralDist(self, p2, border):
        """ Determine lateral distance to a point p2 along profile border.
        """
        path = SegmentedPath()
        p2_project, p2_seg_project = p2.projectOnClosedPath(border)
        project, seg_project = self.projectOnClosedPath(border)
        path.extend([project, p2_project])
        if p2_seg_project < seg_project:
            path.reverse()
        for n in range(min(p2_seg_project, seg_project) + 1,
                       max(p2_seg_project, seg_project)):
            path.insert(len(path)-1, border[n])
        L = path.length()
        return min(L, border.perimeter() - L)

class ProfileBorderData(SegmentedPath):
    def __init__(self, pointlist=[]):
        SegmentedPath.__init__(self, pointlist)
        self.holeli = []

        
    def addHole(self, pointlist=[]):
        self.holeli.append(pointlist)
    
    def area(self):
        """ Determine area of profile, excluding holes
        """
        return SegmentedPath.area(self) - sum([h.area() for h in self.holeli])
        
    def contains(self, p):
        """  Determine whether point p is inside profile border, 
             excluding holes
        """
        if not p:
            return None
        return p.isWithinProfile(self)        
    
class ClusterData(list):
    def __init__(self, pli=[]):
        try:
            self.extend([Particle(p.x, p.y) for p in pli])
        except (AttributeError, IndexError):
            raise TypeError, 'not a particle list'
        self.convexHull = SegmentedPath()

    def lateralDist(self, c2, border):
        """ Determine lateral distance to a cluster c2 along profile border.
        """
        path = SegmentedPath()
        c2_project, c2_seg_project = c2.convexHull.centroid().projectOnClosedPath(border)
        project, seg_project = self.convexHull.centroid().projectOnClosedPath(border)
        path.extend([project, c2_project])
        if c2_seg_project < seg_project:
            path.reverse()
        for n in range(min(c2_seg_project, seg_project) + 1,
                       max(c2_seg_project, seg_project)):
            path.insert(len(path)-1, border[n])
        L = path.length()
        return min(L, border.perimeter() - L)


class ProfileData:
    def __init__(self, inputfn, opt):
        self.inputfn = inputfn
        self.spli = []
        self.lpli = []
        self.randomli = [] 
        self.warnflag = False
        self.errflag = False

    def process(self, opt):
        """ Parse profile data from a file and determine distances
        """
        try:
            self.parse(opt)
            self.checkPaths()
            sys.stdout.write("Processing profile...\n")
            self.posloc = self.path.centroid()
            self.negloc = Point(None, None)
            if not self.posloc.isWithinPolygon(self.path):
                self.negloc = self.posloc
                self.posloc = Point(None, None)
            for p in self.spli:                
                p.determineStuff(self.spli, self, opt)
            self.spli = [p for p in self.spli if not p.skipped]
            for p in self.lpli:                
                p.determineStuff(self.lpli, self, opt)
            self.lpli = [p for p in self.lpli if not p.skipped]
            for r in self.randomli:
                r.determineStuff(self.randomli, self, opt)
            self.randomli = [r for r in self.randomli if not r.skipped]
            self.getInterDistlis(opt)
            self.getClusters(opt)
            self.getMonteCarlo(opt)
            if opt.stop_requested:
                return
            if opt.outputs['individual profiles']:
                self.saveResults(opt)
            sys.stdout.write("Done.\n")
        except ProfileError, (self, msg):
            sys.stdout.write("Error: %s\n" % msg)
            self.errflag = True

    def getMonteCarlo(self, opt):
        if opt.monte_carlo_particle_type == "none":
            return
        sys.stdout.write("Running Monte Carlo simulations:\n")
        self.s_mcli = []
        self.l_mcli = []
        if opt.monte_carlo_particle_type in ("both", "small particles"):
            sys.stdout.write("Simulating small particles...\n")
            self.s_mcli = self.runMonteCarlo("small", opt)
            sys.stdout.write("\n")
        if opt.monte_carlo_particle_type in ("both", "large particles"):
            sys.stdout.write("Simulating large particles...\n")
            self.l_mcli = self.runMonteCarlo("large", opt)
            sys.stdout.write("\n")
            
    def getClusters(self, opt):
        if not (opt.determine_clusters and opt.within_cluster_dist > 0):
            return
        sys.stdout.write("Determining clusters...\n")
        self.spClusters = self.determineClusters(self.spli, opt)
        self.processClusters(self.spClusters, opt)
        self.lpClusters = self.determineClusters(self.lpli, opt)
        self.processClusters(self.lpClusters, opt)            
            
    def getInterDistlis(self, opt):
        if not opt.determine_interpoint_dists:
            return
        if not True in [val for key, val in opt.interpoint_relations.items()
                        if "simulated" not in key]:
            return
        sys.stdout.write("Determining interpoint distances...\n")
        if opt.interpoint_relations["small - small"]:
            self.s_distli, \
                self.s_latdistli = self.getSameInterpointDistances(opt,
                                                                   self.spli)

        if opt.interpoint_relations["large - large"]:
            self.l_distli, \
                self.l_latdistli = self.getSameInterpointDistances(opt,
                                                                   self.lpli)
        if opt.interpoint_relations["small - large"]:
            self.sl_distli, \
                self.sl_latdistli = self.getInterpointDistances2(opt,
                                                                 self.spli,
                                                                 self.lpli)
        if opt.interpoint_relations["large - small"]:
            self.ls_distli,\
                self.ls_latdistli = self.getInterpointDistances2(opt,
                                                                 self.lpli,
                                                                 self.spli)
        if opt.interpoint_relations["random - small"]:
            self.rs_distli, \
                self.rs_latdistli = self.getInterpointDistances2(opt,
                                                                 self.randomli,
                                                                 self.spli)
        if opt.interpoint_relations["random - large"]:
            self.rl_distli, \
                self.rl_latdistli = self.getInterpointDistances2(opt,
                                                                 self.randomli,
                                                                 self.lpli)

    def getSameInterpointDistances(self, opt, pointli):
        dli = []
        latdli = []
        for i in range(0, len(pointli)):
            if opt.stop_requested:
                return [], []
            if opt.interpoint_dist_mode == 'all':
                for j in range(i + 1, len(pointli)):
                    if opt.interpoint_shortest_dist:
                        dli.append(pointli[i].dist(pointli[j]))
                    if opt.interpoint_lateral_dist:
                        latdli.append(pointli[i].lateralDist(pointli[j],
                                                             self.path))
            elif opt.interpoint_dist_mode == 'nearest neighbour':
                if opt.interpoint_shortest_dist:
                    dli.append(pointli[i].determineNearestNeighbour(pointli, self,
                                                                    opt))
                if opt.interpoint_lateral_dist:
                    latdli.append(pointli[i].determineNearestLateralNeighbour(
                                                            pointli, self, opt))
        dli = [d for d in dli if d != None]
        latdli = [d for d in latdli if d != None]
        return dli, latdli


    def getInterpointDistances2(self, opt, pointli, pointli2=[]):
        dli = []
        latdli = []
        for i, p in enumerate(pointli):
            if opt.stop_requested:
                return [], []            
            if opt.interpoint_dist_mode == 'all':
                for p2 in pointli2:
                    if opt.interpoint_shortest_dist:
                        dli.append(p.dist(p2))
                    if opt.interpoint_lateral_dist:
                        latdli.append(p.lateralDist(p2, self.path))
            elif opt.interpoint_dist_mode == 'nearest neighbour':
                if opt.interpoint_shortest_dist:
                    dli.append(p.determineNearestNeighbour(pointli2, self, opt))
                if opt.interpoint_lateral_dist:
                    latdli.append(p.determineNearestLateralNeighbour(pointli2,
                                                                    self, opt))
        dli = [d for d in dli if d != None]
        latdli = [d for d in latdli if d != None]
        return dli, latdli


    def runMonteCarlo(self, ptype, opt):

        def isValid(p):
            if p in mcli[n]["pli"]:
                return False
            if p.isWithinProfile(self.path):
                return True
            if p.isWithinHole(self.path):
                return False
            d = p.perpendDistClosedPath(self.path)
            if d is None:
                return False
            if (opt.monte_carlo_simulation_window == "profile" and not
                opt.monte_carlo_strict_location and
                (abs(d) < toPixelUnits(opt.spatial_resolution,
                                       self.pixelwidth))):
                    return True
            if (opt.monte_carlo_simulation_window == "profile + shell" and
                    abs(d) <= toPixelUnits(opt.shell_width,
                                           self.pixelwidth)):
                return True
            return False

        if ptype == "small":
            if opt.monte_carlo_simulation_window == "profile + shell":
                spli = self.spli    # Points outside the shell have already
                lpli = self.lpli    # been skipped.
            else:                   # If window == "profile"
                spli = [p for p in self.spli if p.isWithinProfile]
                lpli = [p for p in self.lpli if p.isWithinProfile]
            numpoints = len(spli)
        elif ptype == "large":
            if opt.monte_carlo_simulation_window == "profile + shell":
                lpli = self.lpli     # Points outside the shell have already
                spli = self.spli     # been skipped.
            else:                    # If window == "profile"
                lpli = [p for p in self.lpli if p.isWithinProfile]
                spli = [p for p in self.spli if p.isWithinProfile]
            numpoints = len(lpli)
        else:
            return []
        box = self.path.boundingBox()        
        border = 0
        if opt.monte_carlo_simulation_window == "profile + shell":
            border = toPixelUnits(opt.shell_width, self.pixelwidth)
        mcli = []
        for n in range(0, opt.monte_carlo_runs):
            if opt.stop_requested:
                return []            
            dotProgress(n)            
            mcli.append({"pli":[],
                         "same": {"dist": [], "latdist": []},
                         "sim-small": {"dist": [], "latdist": []},
                         "sim-large": {"dist": [], "latdist": []},
                         "small-sim": {"dist": [], "latdist": []},
                         "large-sim": {"dist": [], "latdist": []},
                         "clusterli": []})
            for i in range(0, numpoints):
                while True:
                    x = random.randint(int(box[0].x - border),
                                       int(box[1].x + border) + 1)
                    y = random.randint(int(box[0].y - border),
                                       int(box[2].y + border) + 1)
                    p = Particle(x, y)
                    # escape the while loop when a valid simulated point is found
                    if isValid(p):
                        break
                mcli[n]["pli"].append(p)
            for p in mcli[n]["pli"]:
                p.determineStuff(mcli[n]["pli"], self, opt)
            if opt.determine_interpoint_dists:
                if (opt.interpoint_relations["simulated %s - simulated %s"
                    % (ptype, ptype)]):
                        distlis = self.getSameInterpointDistances(opt,
                                                                mcli[n]["pli"])
                        mcli[n]["same"]["dist"].append(distlis[0])
                        mcli[n]["same"]["latdist"].append(distlis[1])
                if ptype == "small":
                    ptype2 = "large"
                    pli2 = lpli
                else:
                    ptype2 = "small"
                    pli2 = spli
                if (opt.interpoint_relations["simulated %s - %s"
                                            % (ptype, ptype2)]):
                    distlis = self.getInterpointDistances2(opt, mcli[n]["pli"],
                                                           pli2)
                    mcli[n]["sim-%s" % ptype2]["dist"].append(distlis[0])
                    mcli[n]["sim-%s" % ptype2]["latdist"].append(distlis[1])
                if (opt.interpoint_relations["%s - simulated %s"
                                            % (ptype2, ptype)]):
                    distlis = self.getInterpointDistances2(opt, pli2,
                                                           mcli[n]["pli"])
                    mcli[n]["%s-sim" % ptype2]["dist"].append(distlis[0])
                    mcli[n]["%s-sim" % ptype2]["latdist"].append(distlis[1])
        if opt.determine_clusters:
            for n, li in enumerate(mcli):
                mcli[n]["clusterli"] = self.determineClusters(li["pli"], opt)
                self.processClusters(mcli[n]["clusterli"], opt)
        return mcli

    def processClusters(self, clusterli, opt):
        for c in clusterli:
            if opt.stop_requested:
                return                
            c.convexHull = convexHullGraham(c)
            c.distToPath = c.convexHull.centroid(
                                            ).perpendDistClosedPath(self.path)
        for c in clusterli:
            if opt.stop_requested:
                return            
            c.nearestCluster = ClusterData()
            if len(clusterli) == 1:
                c.distToNearestCluster = -1
                return
            c.distToNearestCluster = sys.maxint
            for c2 in clusterli:
                if c2 != c:
                    d = c.lateralDist(c2, self.path)
                    if  d < c.distToNearestCluster:
                        c.distToNearestCluster = d
                        c.nearestCluster = c2


    def determineClusters(self, pointli, opt):
        """ Partition pointli into clusters; each cluster contains all points
            that are less than opt.within_cluster_dist from at least one
            other point in the cluster
        """
        clusterli = []
        for p1 in pointli:
            if opt.stop_requested:
                return []            
            if p1.cluster:
                continue
            for p2 in pointli:
                if p1 != p2 and p1.dist(p2) <= toPixelUnits(
                                                    opt.within_cluster_dist,
                                                    self.pixelwidth):
                    if p2.cluster != None:
                        p1.cluster = p2.cluster
                        clusterli[p1.cluster].append(p1)
                        break
            else:
                p1.cluster = len(clusterli)
                clusterli.append(ClusterData([p1]))
        return clusterli


    def parse(self, opt):
        """ Parse profile data from input file 
        """
        sys.stdout.write("\nParsing '%s':\n" % self.inputfn)
        li = readFile(self.inputfn)
        if not li:
            raise ProfileError(self, "Could not open input file")
        while li:
            str = li.pop(0).replace("\n","").strip()
            if str.split(" ")[0].upper() == "IMAGE":
                self.src_img = str.split(" ")[1]
            elif str.split(" ")[0].upper() == "PROFILE_ID":
                try:
                    self.ID = str.split(" ")[1]
                except IndexError:
                    self.ID = ''
            elif str.split(" ")[0].upper() == "COMMENT":
                try:
                    self.comment = str.split(" ", 1)[1]
                except IndexError:
                    self.comment = ''
            elif str.split(" ")[0].upper() == "PROFILE_TYPE":
                try:
                    self.profile_type = str.split(" ", 1)[1]
                except IndexError:
                    self.profile_type = ''                    
            elif str.split(" ")[0].upper() == "PIXELWIDTH":
                try: 
                    self.pixelwidth = float(str.split(" ")[1])
                    self.metric_unit = str.split(" ")[2]
                except IndexError, ValueError:
                    raise ProfileError(self, 
                                       "PIXELWIDTH is not a valid number")
            elif str.upper() == "PROFILE_BORDER":
                self.path = ProfileBorderData(self.__getCoords(li, "path"))
            elif str.upper() == "PROFILE_HOLE":
                self.path.addHole(SegmentedPath(self.__getCoords(li, 
                                                     "hole"))) 
            elif str.upper() == "SMALL_PARTICLES":
                pli = self.__getCoords(li, "small particle")
                for p in pli: 
                    self.spli.append(Particle(p.x, p.y, type="gold"))
            elif str.upper() == "LARGE_PARTICLES":
                pli = self.__getCoords(li, "large particle")
                for p in pli: 
                    self.lpli.append(Particle(p.x, p.y, type="gold"))   
            elif str.upper() in ("RANDOM", "RANDOM_POINTS"):
                randomli = self.__getCoords(li, "random")
                for r in randomli: 
                    self.randomli.append(Particle(r.x, r.y, type="random"))                    
            elif str[0] != "#":          # unless specifically commented out           
                ProfileWarning(self, "Unrecognized string '" + str + 
                                     "' in input file")
        # Now, let's see if everything was found
        self.checkParsedData(opt)

    def checkParsedData(self, opt):
        """ See if the synapse data was parsed correctly, and print info on the
            parsed data to standard output.            
        """
        try:
            self.src_img
        except AttributeError:
            self.src_img = "N/A"
        sys.stdout.write("  Source image: %s\n" % self.src_img)
        try:
            self.ID
        except AttributeError:
            self.ID = "N/A"
        sys.stdout.write("  Profile ID: %s\n" % self.ID)
        try:
            self.comment
        except AttributeError:
            self.comment = ""
        sys.stdout.write("  Comment: %s\n" % self.comment)
        try:
            self.profile_type
        except AttributeError:
            self.profile_type = ""
        sys.stdout.write("  Profile type: %s\n" % self.comment)        
        try:
            sys.stdout.write("  Pixel width: %.2f " % self.pixelwidth)
        except AttributeError:
            raise ProfileError(self, 
                               "No valid pixel width found in input file")
        try:
            sys.stdout.write("%s\n" % self.metric_unit)
        except AttributeError:
            raise ProfileError(self, "Metric unit not found in input file")
        sys.stdout.write("  Small particles: %d\n" % len(self.spli))
        sys.stdout.write("  Large particles: %d\n" % len(self.lpli))
        try:
            self.path[0]
        except (IndexError, AttributeError):
            raise ProfileError(self, "Profile border not found in input file")
        sys.stdout.write("  Profile holes: %d\n" % len(self.path.holeli))            
        try:
            self.randomli[0]
        except (IndexError, AttributeError):
            try:
                if opt.use_random:
                    raise ProfileError(self, "Random points expected but not"
                                             " found")
            except AttributeError:
                opt.use_random = False
        else:
            try:
                if not opt.use_random:
                    raise ProfileError(self, "Random points found here but not"
                                             " in previous input files")
            except AttributeError:
                opt.use_random = True
            sys.stdout.write("  Random points specified.\n")

    def checkPaths(self):
        def checkPath(path, s):
            for n1 in range(0, len(path)-3):
                for n2 in range(0, len(path)-1):
                    if n1 not in (n2, n2+1) and n1+1 not in (n2, n2+1):
                        if segmentIntersection(path[n1], path[n1+1],
                                               path[n2], path[n2+1]):
                            raise ProfileError(self,
                                            "%s invalid (crosses itself)" % s)
            return True

        checkPath(self.path, "Profile border")
        for path in self.path.holeli:
            checkPath(path, "Hole")
        for n, h in enumerate(self.path.holeli):
            if not h.isSimplePolygon():
                raise ProfileError(self,
                                   "Profile hole %d is not a simple polygon"
                                    % (n+1))
            if not h.isWithinPolygon(self.path):
                raise ProfileError(self,
                                   "Profile hole %d is not (completely) "
                                   "within profile" % (n+1))
            for n2, h2 in enumerate(self.path.holeli[n+1:]):
                if h.overlapsPolygon(h2):
                    raise ProfileError(self,
                                       "Profile hole %d overlaps with hole %d "
                                       % (n+1, n+n2+2))
        sys.stdout.write("  Paths are ok.\n")


    def __getCoords(self, strli, coordType=""):
        """ Pop point coordinates from list strli.
            When an element of strli is not a valid point,
            a warning is issued.
        """
        pointli = []
        s = strli.pop(0).replace("\n","").replace(" ","").strip()
        while s != "END":
            try:
                p = Point(float(s.split(",")[0]), float(s.split(",")[1]))
                if pointli and (p == pointli[-1] or 
                                (coordType in ('particle', 'random')
                                 and p in pointli)):
                    sys.stdout.write("Duplicate %s coordinates %s: skipping "
                                     "2nd instance\n" % (coordType, p))
                else:
                    pointli.append(p)                    
            except ValueError:
                if s[0] != "#":
                    ProfileWarning(self, "'%s' not valid %s coordinates" 
                                   % (s, coordType))
                else:
                    pass 
            s = strli.pop(0).replace("\n","").strip()
        # For some reason, sometimes the endnodes have the same coordinates;
        # in that case, delete the last endnode to avoid division by zero
        if (len(pointli) > 1) and (pointli[0] == pointli[-1]):
            del pointli[-1]
        return pointli        

    def saveResults(self, opt):
        """ Output results from a single synapse to file
        """ 
        
        def m(x):
            try:
                return toMetricUnits(x, self.pixelwidth)
            except ZeroDivisionError:
                return None
                
        def m2(x): 
            try:
                return toMetricUnits(x, self.pixelwidth**2) # for area units...
            except ZeroDivisionError:
                return None
       
        def fwrite(*args):
            f.writerow(args)
                    
        
        try:
            self.outputfn = os.path.join(opt.outputDir,
                                         opt.output_filename_suffix +
                                         os.path.basename(self.inputfn)
                                         + opt.output_filename_ext)
            if (os.path.exists(self.outputfn) and
                opt.action_if_output_file_exists == 'enumerate'):
                    self.outputfn = enumFilename(self.outputfn, 2)
            sys.stdout.write("Writing to '%s'...\n" % self.outputfn)
            if opt.outputFileFormat == "csv":
                f = csv.writer(file(self.outputfn, "w"), **opt.csv_format)
            elif opt.outputFileFormat == 'excel':
                import xls
                f = xls.writer(self.outputfn)
            fwrite("Table 1. Profile-centric data")
            fwrite("Source image:", self.src_img)
            fwrite("Profile ID:", self.ID)
            fwrite("Comment:", self.comment)
            fwrite("Profile type:", self.profile_type) 
            fwrite("Pixel width:", tostr(float(self.pixelwidth), 2), 
                                   self.metric_unit)
            fwrite("Number of small particles (total):", len(self.pli))            
            fwrite("Number of large particles (total):", len(self.pli))
            fwrite("Number of random points (total):", len(self.randomli))
            fwrite("Table 2. Small particle-centric data")                                   
            columnheadings = ["Particle number (as appearing in input file)",
                              "Particle coordinates (in pixels)"]
            fwrite(*columnheadings) 
            f.writerows([[n+1,
                          str(p)] 
                          for n, p in enumerate(self.spli)])     
            fwrite("Table 3. Large particle-centric data")                                   
            columnheadings = ["Particle number (as appearing in input file)",
                              "Particle coordinates (in pixels)"]
            fwrite(*columnheadings) 
            f.writerows([[n+1,
                          str(p)] 
                          for n, p in enumerate(self.lpli)])                               
            fwrite("Table 4. Random point-centric data")
            columnheadings = ["Random point number (as appearing in input file)",
                              "Random point coordinates (in pixels)"]
            fwrite(*columnheadings) 
            f.writerows([[n+1,
                          str(r)] 
                          for n, r in enumerate(self.randomli)])                               
            f.close()
        except IOError:
            raise ProfileError(self, "Unable to write to file '%s'" 
                               % self.outputfn)
        sys.stdout.write("Done.\n")
        return 1


# end of class Profile        
        
class OptionData:
    def __init__(self):
        self.input_file_list = []
        self.spatial_resolution = 25
        self.shell_width = 0   # Skip points farther than this from profile
        self.outputs = {'profile summary': True, 'point summary': True,
                        'random summary': True, 'session summary': True,
                        'individual profiles': False}
        self.output_file_format = 'excel'
        self.output_filename_ext = '.xls'
        self.input_filename_ext = '.pd'
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

def ProfileWarning(profile, msg):
    """ Issue a warning
    """
    sys.stdout.write("Warning: %s.\n" % msg)
    profile.warnflag = True      

def ProfileMessage(profile, msg):
    """ Show a message
    """
    sys.stdout.write("%s.\n" % msg)

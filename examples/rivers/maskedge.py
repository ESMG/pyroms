
import numpy as np
import netCDF4
import sys
#import Ngl

# This thing was originally written in C++ before there was a standard
# (1994 or 1995). Unfortunately, it needed to be recompiled for each
# domain so here I attempt a rewrite.
#
# There was even once a version that used the graphics from the HLU
# library that was a step on the road to NCL. Hmm - plots... include as
# many of the Ngl commands as you like.
#                                 - Kate Hedstrom (2015)

# Read the rho-point land mask from a NetCDF grid file and produce
# the i,j pairs describing the land/sea interface. This was
# originally for the elliptic solver in SPEM, but is now for river
# sources in ROMS. It can also be used just to find out if you left
# any lakes in your domain when editing the land mask.

def readmask(ncfile):
    nc = netCDF4.Dataset(ncfile, 'r')

    mask2 = nc.variables['mask_rho'][:]
    imask = mask2.astype(int)
    nc.close()

    # Stop the flood-fill at the edges
    imask[0,:] = -1
    imask[:,0] = -1
    imask[-1,:] = -1
    imask[:,-1] = -1
    return imask

def flood_fill_water(imask, i, j, ii):
    """
    Find all water touching point (i,j) and set to ii.
    The recursive version works for teeny, tiny grids only, so we need
    another way.
    """
    Mp, Lp = imask.shape
    Lm = Lp-2
    Mm = Mp-2
    jj = imask[j,i]

    llist = []
    llist.append((j,i))
    while len(llist) > 0:
        (j,i) = llist.pop()
        imask[j,i] = ii
        if ( imask[j,i-1] == jj and i > 1 ):
            llist.append((j, i-1))
        if ( imask[j-1,i] == jj and j > 1 ):
            llist.append((j-1, i))
        if ( imask[j,i+1] == jj and i < Lm ):
            llist.append((j, i+1))
        if ( imask[j+1,i] == jj and j < Mm ):
            llist.append((j+1, i))

def flood_fill_land(imask, i, j, ii):
    """
    Find all land touching point (i,j) and set to ii.
    The recursive version works for teeny, tiny grids only, so we need
    another way.
    """
    Mp, Lp = imask.shape
    Lm = Lp-2
    Mm = Mp-2
    jj = imask[j,i]

    llist = []
    llist.append((j,i))
    while len(llist) > 0:
        (j,i) = llist.pop()
        imask[j,i] = ii
        if ( imask[j,i-1] == jj and i > 1 ):
            llist.append((j, i-1))
        if ( imask[j-1,i] == jj and j > 1 ):
            llist.append((j-1, i))
        if ( imask[j,i+1] == jj and i < Lm ):
            llist.append((j, i+1))
        if ( imask[j+1,i] == jj and j < Mm ):
            llist.append((j+1, i))
# now do the diagonals
        if ( imask[j-1,i-1] == jj and i > 1 ):
            llist.append((j-1, i-1))
        if ( imask[j-1,i+1] == jj and j > 1 ):
            llist.append((j-1, i+1))
        if ( imask[j+1,i+1] == jj and i < Lm ):
            llist.append((j+1, i+1))
        if ( imask[j+1,i-1] == jj and j < Mm ):
            llist.append((j+1, i-1))

def set_values(imask, k, val):
    """ Set all k values to val"""
    Mp, Lp = imask.shape
    for j in range(1,Mp-1):
        for i in range(1,Lp-1):
            if (imask[j,i] == k): imask[j,i] = val

def warning(*objs):
    print("STDERR: ", *objs, file=sys.stderr)

def color_water(imask):
    """ Set each body of water to a different value, return number of
    bodies of water."""
    count = 2

    Mp, Lp = imask.shape
    for j in range(1,Mp-1):
        for i in range(1,Lp-1):
            if (imask[j,i] == 1):
                flood_fill_water(imask, i, j, count)
                warning("New body!", i, j)
                count += 1

    warning("There are", count-2, " bodies of water.")
    return count-2

def peninsula(imask, plist, i, j, dir, iwat, iland):
    """
    Process a found peninsula - append to a linked list of its points
    and then set its mask points to water.
    """
    Mp, Lp = imask.shape
    plist.append((-3, -3))
    p = (i,j)

    # set up new linked list of Points
#    assert(isEdge(p));
    plist.append(p)

    if (dir == 'east'):
        seed = (i,j)
    elif (dir == 'west'):
        seed = (i-1,j-1)
    elif (dir == 'south'):
        seed = (i,j-1)
    elif (dir == 'north'):
        seed = (i-1,j)
#    warning("Peninsula at", seed, dir)

    # Trace the edge of the peninsula, keeping track of the psi
    # points in a linked list.  The seed point is kept for the
    # flood_fill, setting the land to match the current water
    # value.
    while True:
        if (dir == 'east'):
            i += 1
            p = (i,j)
            if ((imask[j-1,i] == iland) and (imask[j-1,i-1] == iwat)):
                dir = 'south'
            elif ((imask[j,i] == iland) and (imask[j-1,i] == iwat)):
                dir = 'east'
            elif ((imask[j,i-1] == iland) and (imask[j,i] == iwat)):
                dir = 'north'
            elif (i==1 or j==1 or i==Lp-1 or j==Mp-1):
                break
            else:
                warning("Problem in peninsula at ", i, j)
                exit(1)
        elif (dir == 'north'):
            j += 1
            p = (i,j)
            if ((imask[j,i] == iland) and (imask[j-1,i] == iwat)):
                dir = 'east'
            elif ((imask[j,i-1] == iland) and (imask[j,i] == iwat)):
                dir = 'north'
            elif ((imask[j-1,i-1] == iland) and (imask[j,i-1] == iwat)):
                dir = 'west'
            elif (i==1 or j==1 or i==Lp-1 or j==Mp-1):
                break
            else:
                warning("Problem in peninsula at ", i, j)
                exit(1)
        elif (dir == 'west'):
            i -= 1
            p = (i,j)
            if ((imask[j,i-1] == iland) and (imask[j,i] == iwat)):
                dir = 'north'
            elif ((imask[j-1,i-1] == iland) and (imask[j,i-1] == iwat)):
                dir = 'west'
            elif ((imask[j-1,i] == iland) and (imask[j-1,i-1] == iwat)):
                dir = 'south'
            elif (i==1 or j==1 or i==Lp-1 or j==Mp-1):
                break
            else:
                warning("Problem in peninsula at ", i, j)
                exit(1)
        elif (dir == 'south'):
            j -= 1
            p = (i,j)
            if ((imask[j-1,i-1] == iland) and (imask[j,i-1] == iwat)):
                dir = 'west'
            elif ((imask[j-1,i] == iland) and (imask[j-1,i-1] == iwat)):
                dir = 'south'
            elif ((imask[j,i] == iland) and (imask[j-1,i] == iwat)):
                dir = 'east'
            elif (i==1 or j==1 or i==Lp-1 or j==Mp-1):
                break
            else:
                warning("Problem in peninsula at ", i, j)
                exit(1)

        plist.append(p)
    plist.append(p)

    # time to flood_fill to show we are done with this peninsula
    i = seed[0];
    j = seed[1];
    flood_fill_land(imask, i, j, iwat);

def island(imask, ilist, i, j, dir, iwat, iland):
    """
    Process a found island - create a linked list of its points
    and then set its mask points to water.
    """
    ilist.append((-1, -1))
    p = (i,j)
    pstart = p

    # add to list of i,j Points
#    assert(!p.isEdge());
    ilist.append(p)
    # Now we change p so the exit condition doesn't happen yet...
    p = (0,0)

    if (dir == 'east'):
        seed = (i,j)
    elif (dir == 'west'):
        seed = (i-1,j-1)
    elif (dir == 'south'):
        seed = (i,j-1)
    elif (dir == 'north'):
        seed = (i-1,j)
#    warning("Island at", seed, dir)

    # Trace the edge of the island, keeping track of the psi
    # points in a linked list.  Also keep track of the east and west
    # edge segments so that we can later change the peninsula mask
    # points to water.
    while True:
        if (p == pstart):
            break
        if (dir == 'east'):
            i += 1
            p = (i,j)
            if ((imask[j-1,i] == iland) and (imask[j-1,i-1] == iwat)):
                dir = 'south'
            elif ((imask[j,i] == iland) and (imask[j-1,i] == iwat)):
                dir = 'east'
            elif ((imask[j,i-1] == iland) and (imask[j,i] == iwat)):
                dir = 'north'
            else:
                warning("Problem in island at ", i, j)
                exit(1)
        elif (dir == 'north'):
            j += 1
            p = (i,j)
            if ((imask[j,i] == iland) and (imask[j-1,i] == iwat)):
                dir = 'east'
            elif ((imask[j,i-1] == iland) and (imask[j,i] == iwat)):
                dir = 'north'
            elif ((imask[j-1,i-1] == iland) and (imask[j,i-1] == iwat)):
                dir = 'west'
            else:
                warning("Problem in island at ", i, j)
                exit(1)
        elif (dir == 'west'):
            i -= 1
            p = (i,j)
            if ((imask[j,i-1] == iland) and (imask[j,i] == iwat)):
                dir = 'north'
            elif ((imask[j-1,i-1] == iland) and (imask[j,i-1] == iwat)):
                dir = 'west'
            elif ((imask[j-1,i] == iland) and (imask[j-1,i-1] == iwat)):
                dir = 'south'
            else:
                warning("Problem in island at ", i, j)
                exit(1)
        elif (dir == 'south'):
            j -= 1
            p = (i,j)
            if ((imask[j-1,i-1] == iland) and (imask[j,i-1] == iwat)):
                dir = 'west'
            elif ((imask[j-1,i] == iland) and (imask[j-1,i-1] == iwat)):
                dir = 'south'
            elif ((imask[j,i] == iland) and (imask[j-1,i] == iwat)):
                dir = 'east'
            else:
                warning("Problem in island at ", i, j)
                exit(1)

        ilist.append(p)

  # time to flood_fill to show we are done with this island
    i = seed[0];
    j = seed[1];
    flood_fill_land(imask, i, j, iwat);

def edges(imask, plist, iwat, iland):
    """
    Search for peninsulas by going around the boundary.
    """
    Mp, Lp = imask.shape
    for j in range(1,Mp-2):
        if ((imask[j,1] == iwat) and (imask[j+1,1] == iland)):
            peninsula(imask, plist, 1, j+1, 'east', iwat, iland)
    for i in range(1,Lp-2):
        if ((imask[Mp-2,i] == iwat) and (imask[Mp-2,i+1] == iland)):
            peninsula(imask, plist, i+1, Mp-1, 'south', iwat, iland)
    for j in list(reversed(list(range(2,Mp-1)))):
        if ((imask[j,Lp-2] == iwat) and (imask[j-1,Lp-2] == iland)):
            peninsula(imask, plist, Lp-1, j, 'west', iwat, iland)
    for i in list(reversed(list(range(2,Lp-1)))):
        if ((imask[1,i] == iwat) and (imask[1,i-1] == iland)):
            peninsula(imask, plist, i, 1, 'north', iwat, iland)

def interior(imask, ilist, iwat, iland):
    """
    Search for islands by scanning the interior.
    """
    Mp, Lp = imask.shape
    for i in range(2,Lp-2):
        for j in range(2,Mp-2):
            if ((imask[j,i] == iwat) and (imask[j+1,i] == iland)):
                island(imask, ilist, i, j+1, 'east', iwat, iland)

def main():
    ncfile = sys.argv[1]
    imask = readmask(ncfile)
    lpoints = []
    ipoints = []

#    wks_type = "X11"
#    wks = Ngl.open_wks(wks_type,"maskedge")
#    res = Ngl.Resources()
#    res.cnFillOn = True

#    Ngl.contour(wks,imask,res)
    count = color_water(imask)
#    Ngl.contour(wks,imask,res)
    iland = 0
    for iwat in range(2, count+2):
        edges(imask, lpoints, iwat, iland)
#        Ngl.contour(wks,imask,res)
        interior(imask, ipoints, iwat, iland)
#        Ngl.contour(wks,imask,res)
        set_values(imask, iland, iwat)
#        Ngl.contour(wks,imask,res)
        iland = iwat

    # Islands first, then peninsulas
    for point in ipoints:
        i,j = point
        print(i, j)
    for point in lpoints:
        i,j = point
        print(i, j)
    print(-10, -10)
#    Ngl.end()

if __name__ == "__main__":
    main()

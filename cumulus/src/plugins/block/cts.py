
# This file is part of the Coriolis Software.
# Copyright (c) Sorbonne Université 2026-2026, All Rights Reserved
#
# +-----------------------------------------------------------------+
# |                   C O R I O L I S                               |
# |      C u m u l u s  -  P y t h o n   T o o l s                  |
# |                                                                 |
# |  Author      :                           Hongda YANG            |
# |  E-mail      :                  laskcl4116@gmail.com            |
# | =============================================================== |
# |  Python      :       "./plugins/block/cts.py"                   |
# +-----------------------------------------------------------------+

"""
Post-placement buffering utilities, building blocks for a clock tree
synthesis performed *after* the placement.

The placement must be completed and EtesianEngine.toHurricane() must have
been called, so the rows are filled with feed cells. Buffers are inserted
by replacing feed cells (see EtesianEngine.insertCellNear()), so the
placement of the other cells is never modified.

Only flat designs are supported for now: the net to buffer must belong
to the cell placed by Etesian (no block instance).
"""

from ...Hurricane       import Box, Point, Net, RoutingPad, Occurrence
from ...helpers         import trace
from ...helpers.io      import ErrorMessage
from ...helpers.overlay import UpdateSession


def getBufferIos ( bufferCell ):
    """
    Return the (input, output) master nets of a buffer cell, that is, its
    only non-supply input and its only output.
    """
    inputs  = []
    outputs = []
    for net in bufferCell.getNets():
        if not net.isExternal(): continue
        if net.isSupply(): continue
        direction = net.getDirection() & Net.Direction.INOUT
        if   direction == Net.Direction.IN:  inputs .append( net )
        elif direction == Net.Direction.OUT: outputs.append( net )
    if len(inputs) != 1 or len(outputs) != 1:
        raise ErrorMessage( 1, [ 'cts.getBufferIos(): "{}" is not a buffer.'.format( bufferCell.getName() )
                               , 'It must have exactly one input and one output '
                                 '(found {} and {}).'.format( len(inputs), len(outputs) ) ] )
    return inputs[0], outputs[0]


def getSinksCenter ( sinkPlugs ):
    """Return the center of gravity of the instances of ``sinkPlugs``."""
    x = 0
    y = 0
    for plug in sinkPlugs:
        center = plug.getInstance().getAbutmentBox().getCenter()
        x += center.getX()
        y += center.getY()
    return Point( x // len(sinkPlugs), y // len(sinkPlugs) )


def placeBufferNear ( etesian, bufferCell, name, position, sliceHeight, maxRows=10 ):
    """
    Insert ``bufferCell`` in the placed rows, as close as possible to
    ``position``. The row of ``position`` is tried first, then the rows
    above and below it, alternately, up to ``maxRows`` rows away. In each
    row, the feed closest to the requested abscissa is used.

    Return the new instance, or None if no room has been found.
    """
    cellAb = etesian.getCell().getAbutmentBox()
    xHint  = position.getX() - bufferCell.getAbutmentBox().getWidth() // 2
    for distance in range( maxRows+1 ):
        for side in ((0,) if distance == 0 else (1,-1)):
            y = position.getY() + side*distance*sliceHeight
            if y < cellAb.getYMin() or y >= cellAb.getYMax(): continue
            area     = Box( cellAb.getXMin(), y - sliceHeight//4
                          , cellAb.getXMax(), y + sliceHeight//4 )
            instance = etesian.insertCellNear( bufferCell, name, area, xHint )
            if instance:
                trace( 550, '\tcts.placeBufferNear(): {} at {}\n'.format( instance, area ))
                return instance
    return None


def insertBuffer ( etesian, net, sinkPlugs, bufferCell, name, sliceHeight, maxRows=10 ):
    """
    Insert, after placement, a buffer driving a subset of the sinks of ``net``.

    Before::

        net --+--> sink 1
              +--> sink 2
              +--> other sinks

    After::

        net --+--> buffer --> <name>_out --+--> sink 1
              |                            +--> sink 2
              +--> other sinks

    The buffer is placed as close as possible to the center of gravity of
    ``sinkPlugs``, in a feed cell (see placeBufferNear()). The output net
    gets the same type as ``net``, so a clock stays a clock.

    If ``net`` has already been flattened (it has RoutingPads), the
    RoutingPads of the moved sinks are moved to the output net, and new
    ones are created for the buffer terminals. Otherwise (the clock nets
    are usually not flattened before routing), nothing is done, the router
    will flatten both nets.

    Return a tuple (bufferInstance, outputNet).
    """
    cell = net.getCell()
    if cell != etesian.getCell():
        raise ErrorMessage( 1, [ 'cts.insertBuffer(): Net "{}" does not belong to the placed cell "{}".' \
                                 .format( net.getName(), etesian.getCell().getName() )
                               , '(hierarchical designs are not supported yet)' ] )
    if not sinkPlugs:
        raise ErrorMessage( 1, 'cts.insertBuffer(): No sink given for net "{}".'.format( net.getName() ))
    for plug in sinkPlugs:
        if plug.getNet() != net:
            raise ErrorMessage( 1, 'cts.insertBuffer(): {} is not connected to net "{}".' \
                                   .format( plug, net.getName() ))
        if plug.getMasterNet().getDirection() & Net.Direction.OUT:
            raise ErrorMessage( 1, 'cts.insertBuffer(): {} is not a sink (it is an output).'.format( plug ))
    outputName = '{}_out'.format( name )
    if cell.getNet( outputName ):
        raise ErrorMessage( 1, 'cts.insertBuffer(): A net named "{}" already exists in "{}".' \
                               .format( outputName, cell.getName() ))
    bufferInput, bufferOutput = getBufferIos( bufferCell )

    with UpdateSession():
        instance = placeBufferNear( etesian, bufferCell, name, getSinksCenter(sinkPlugs)
                                  , sliceHeight, maxRows )
        if instance is None:
            return None, None

        flattened = False
        for rp in net.getRoutingPads():
            flattened = True
            break

        outputNet = Net.create( cell, outputName )
        outputNet.setType( net.getType() )
        inputPlug  = instance.getPlug( bufferInput  )
        outputPlug = instance.getPlug( bufferOutput )
        inputPlug .setNet( net       )
        outputPlug.setNet( outputNet )

        movedIds = set( plug.getId() for plug in sinkPlugs )
        if flattened:
            staleRps = [ rp for rp in net.getRoutingPads()
                         if rp.getPlugOccurrence().getEntity().getId() in movedIds ]
            for rp in staleRps:
                if len( list(rp.getSlaveComponents()) ):
                    raise ErrorMessage( 1, 'cts.insertBuffer(): {} is already routed.'.format( rp ))
                rp.destroy()
        for plug in sinkPlugs:
            plug.setNet( outputNet )
        if flattened:
            for plug in sinkPlugs:
                RoutingPad.create( outputNet, Occurrence(plug), RoutingPad.BiggestArea )
            RoutingPad.create( net      , Occurrence(inputPlug ), RoutingPad.BiggestArea )
            RoutingPad.create( outputNet, Occurrence(outputPlug), RoutingPad.BiggestArea )

    trace( 550, '\tcts.insertBuffer(): {} drives {} sinks of "{}" through "{}".\n' \
                .format( instance, len(sinkPlugs), net.getName(), outputName ))
    return instance, outputNet


def getSinks ( net ):
    """Return the plugs of ``net`` which are not outputs (its sinks)."""
    return [ plug for plug in net.getPlugs()
             if not (plug.getMasterNet().getDirection() & Net.Direction.OUT) ]


def partitionSinks ( sinkPlugs, maxSinks ):
    """
    Split ``sinkPlugs`` in groups of at most ``maxSinks`` geographically
    close plugs, by recursive bisection. At each step, the plugs are sorted
    along the larger side of their bounding box and cut so that both sides
    get a number of groups proportional to their number of plugs. All the
    groups end up with nearly the same size.

    Return a list of lists of plugs.
    """
    groupCount = (len(sinkPlugs) + maxSinks - 1) // maxSinks
    if groupCount <= 1:
        return [ sinkPlugs ]
    centers = {}
    bb      = Box()
    for plug in sinkPlugs:
        center = plug.getInstance().getAbutmentBox().getCenter()
        centers[ plug.getId() ] = center
        bb.merge( center )
    if bb.getWidth() >= bb.getHeight():
        key = lambda plug: (centers[plug.getId()].getX(), centers[plug.getId()].getY())
    else:
        key = lambda plug: (centers[plug.getId()].getY(), centers[plug.getId()].getX())
    ordereds   = sorted( sinkPlugs, key=key )
    leftGroups = groupCount // 2
    cut        = (len(ordereds) * leftGroups) // groupCount
    return partitionSinks( ordereds[:cut], maxSinks ) \
         + partitionSinks( ordereds[cut:], maxSinks )


def buildClockTree ( etesian, net, bufferCell, sliceHeight, maxSinks, maxRows=10 ):
    """
    Build, after placement, a buffer tree on ``net`` so that no net of the
    tree drives more than ``maxSinks`` sinks. The tree is built bottom-up:
    the sinks of ``net`` are partitioned in groups (see partitionSinks()),
    each group gets its own buffer (see insertBuffer()), then the inputs
    of those buffers, which are now the sinks of ``net``, are partitioned
    in turn, until ``net`` itself has at most ``maxSinks`` sinks.

    Buffers are named "cts_<net>_<level>_<index>", level 0 being the one
    closest to the original sinks.

    Return the list of the levels, each level being the list of its
    buffer instances.
    """
    if maxSinks < 2:
        raise ErrorMessage( 1, 'cts.buildClockTree(): maxSinks must be at least 2 (got {}).'.format( maxSinks ))
    levels = []
    while True:
        sinks = getSinks( net )
        if len(sinks) <= maxSinks: break
        level = []
        for group in partitionSinks( sinks, maxSinks ):
            name = 'cts_{}_{}_{}'.format( net.getName(), len(levels), len(level) )
            instance, outputNet = insertBuffer( etesian, net, group, bufferCell, name
                                              , sliceHeight, maxRows )
            if instance is None:
                raise ErrorMessage( 1, [ 'cts.buildClockTree(): No room left to place buffer "{}" of "{}".' \
                                         .format( name, net.getName() )
                                       , 'Searched {} rows around the sinks center, please increase free space.' \
                                         .format( maxRows ) ] )
            level.append( instance )
        levels.append( level )
    return levels

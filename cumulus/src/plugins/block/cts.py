
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

from ...Hurricane       import DbU, Box, Point, Net, RoutingPad, Occurrence
from ...CRL             import AllianceFramework, Catalog
from ...helpers         import trace
from ...helpers.io      import ErrorMessage, WarningMessage
from ...helpers.overlay import UpdateSession
from .nldm              import NldmLibrary


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
    ``position``. The search is done in a window around ``position``,
    starting at two rows wide and doubled until room is found, up to the
    whole width of the cell. In a window, the row of ``position`` is tried
    first, then the rows above and below it, alternately, as long as they
    are inside the window and at most ``maxRows`` rows away. In each row,
    the feed closest to the requested abscissa is used. So a near feed in
    the next row is preferred over a far one in the row of ``position``.

    Return the new instance, or None if no room has been found.
    """
    cellAb = etesian.getCell().getAbutmentBox()
    xHint  = position.getX() - bufferCell.getAbutmentBox().getWidth() // 2
    window = 2*sliceHeight
    while True:
        for distance in range( maxRows+1 ):
            if distance*sliceHeight > window: break
            for side in ((0,) if distance == 0 else (1,-1)):
                y = position.getY() + side*distance*sliceHeight
                if y < cellAb.getYMin() or y >= cellAb.getYMax(): continue
                area     = Box( max( cellAb.getXMin(), xHint-window ), y - sliceHeight//4
                              , min( cellAb.getXMax(), xHint+window ), y + sliceHeight//4 )
                instance = etesian.insertCellNear( bufferCell, name, area, xHint )
                if instance:
                    trace( 550, '\tcts.placeBufferNear(): {} at {}\n'.format( instance, area ))
                    return instance
        if window >= cellAb.getWidth(): return None
        window *= 2


def insertBuffer ( etesian, net, sinkPlugs, bufferCell, name, sliceHeight, maxRows=10, position=None ):
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

    The buffer is placed as close as possible to ``position``, by default
    the center of gravity of ``sinkPlugs``, in feed cells (see
    placeBufferNear()). The output net
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
        if position is None:
            position = getSinksCenter( sinkPlugs )
        instance = placeBufferNear( etesian, bufferCell, name, position, sliceHeight, maxRows )
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


def getTreeLength ( points ):
    """Length of the Manhattan minimum spanning tree of ``points`` (Prim)."""
    if len(points) < 2: return 0
    inTree   = [ False ] * len(points)
    distance = [ None  ] * len(points)
    distance[0] = 0
    length      = 0
    for _ in range( len(points) ):
        i = min( (j for j in range(len(points)) if not inTree[j] and distance[j] is not None)
               , key=lambda j: distance[j] )
        inTree[i] = True
        length   += distance[i]
        for j in range( len(points) ):
            if inTree[j]: continue
            d = abs(points[i].getX()-points[j].getX()) + abs(points[i].getY()-points[j].getY())
            if distance[j] is None or d < distance[j]: distance[j] = d
    return length


class BufferSizer ( object ):
    """
    Choose the buffer of each group of sinks from its Liberty (NLDM) data:
    the smallest of ``buffers`` whose output transition, under its estimated
    load, does not exceed ``maxSlew`` (the input transition is assumed to be
    ``maxSlew`` too), within its max_capacitance. If none fits, the strongest
    one is used.

    The load is the sum of the sinks input capacitances plus the wire,
    estimated as the minimum spanning tree of the sinks and of their center
    of gravity (where the buffer is placed) times ``wireCapacitance``.

    ``maxSlew`` and the capacitances are in the units of the Liberty file
    (usually ns and pF), ``wireCapacitance`` being per micron. ``liberty``
    defaults to the one given to Yosys by the PDK setup.

    If ``wireResistance`` (ohm per micron) is given, the sizer also balances
    the delays, in the spirit of the Deferred-Merge Embedding: the buffer of
    a group is not put at the center of gravity of its sinks but where the
    sums of the wire delay to each sink and of the delay already below that
    sink (zero for a flip-flop, the estimated latency of a lower level
    buffer) are the closest (see getPosition()).

    An instance can be given as the ``bufferCells`` of buildClockTree().
    """

    GridSize = 25

    def __init__ ( self, buffers, maxSlew, wireCapacitance, liberty=None, wireResistance=None ):
        if liberty is None:
            from ...designflow.yosys import Yosys
            liberty = Yosys._liberty
        if liberty is None:
            raise ErrorMessage( 1, 'cts.BufferSizer(): No Liberty file given and none set by the PDK.' )
        self.library         = NldmLibrary( str(liberty) )
        self.maxSlew         = maxSlew
        self.wireCapacitance = wireCapacitance
        self.wireResistance  = wireResistance
        self.latencies       = {}
        self.missings        = set()
        framework            = AllianceFramework.get()
        self.buffers         = []
        for name in buffers:
            cell   = framework.getCell( name, Catalog.State.Views )
            timing = self.library.getCell( name )
            if cell is None or timing is None:
                raise ErrorMessage( 1, 'cts.BufferSizer(): Buffer "{}" not found in the {}.' \
                                       .format( name, 'library' if cell is None else 'Liberty file' ))
            self.buffers.append( (cell, timing) )
        self.buffers.sort( key=lambda item: sum(item[1].pins.values()) )

    def getPinCapacitance ( self, plug ):
        name   = str( plug.getInstance().getMasterCell().getName() )
        timing = self.library.getCell( name )
        if timing is None:
            if name not in self.missings:
                print( WarningMessage( 'cts.BufferSizer(): No Liberty data for "{}", '
                                       'its input capacitance is ignored.'.format(name) ))
                self.missings.add( name )
            return 0.0
        return timing.pins.get( str(plug.getMasterNet().getName()), 0.0 )

    def getLoad ( self, sinkPlugs, position=None ):
        if position is None: position = getSinksCenter( sinkPlugs )
        points = [ position ] + [ plug.getInstance().getAbutmentBox().getCenter() for plug in sinkPlugs ]
        wire   = DbU.toPhysical( getTreeLength(points), DbU.UnitPowerMicro ) * self.wireCapacitance
        return wire + sum( self.getPinCapacitance(plug) for plug in sinkPlugs )

    def getArrivals ( self, sinkPlugs, position ):
        """
        Estimated delay from ``position`` to the leaves below each sink: the
        Elmore delay of a direct wire to the sink plus the latency already
        recorded for its instance (see setLatency()).
        """
        arrivals = []
        for plug in sinkPlugs:
            center = plug.getInstance().getAbutmentBox().getCenter()
            length = DbU.toPhysical( abs(center.getX()-position.getX()) + abs(center.getY()-position.getY())
                                   , DbU.UnitPowerMicro )
            wire   = self.wireResistance * length \
                   * (self.wireCapacitance*length/2 + self.getPinCapacitance(plug)) * 1e-3
            arrivals.append( wire + self.latencies.get( plug.getInstance().getId(), 0.0 ))
        return arrivals

    def getPosition ( self, sinkPlugs ):
        """
        The point, on a GridSize x GridSize grid over the bounding box of the
        sinks (plus their center of gravity), where the spread of the
        arrivals (see getArrivals()) is the smallest, ties broken by the
        total wire length. Without ``wireResistance``, the center of gravity.
        """
        center = getSinksCenter( sinkPlugs )
        if self.wireResistance is None: return center
        bb = Box()
        for plug in sinkPlugs: bb.merge( plug.getInstance().getAbutmentBox().getCenter() )
        candidates = [ center ]
        for i in range( self.GridSize ):
            for j in range( self.GridSize ):
                candidates.append( Point( bb.getXMin() + (bb.getWidth () * i) // (self.GridSize-1)
                                        , bb.getYMin() + (bb.getHeight() * j) // (self.GridSize-1) ))
        def cost ( point ):
            arrivals = self.getArrivals( sinkPlugs, point )
            length   = sum( abs(p.getInstance().getAbutmentBox().getCenter().getX()-point.getX())
                          + abs(p.getInstance().getAbutmentBox().getCenter().getY()-point.getY())
                            for p in sinkPlugs )
            return ( round(max(arrivals) - min(arrivals), 4), length )
        return min( candidates, key=cost )

    def setLatency ( self, instance, sinkPlugs ):
        """
        Record the estimated latency of the freshly inserted buffer ``instance``
        driving ``sinkPlugs``: its cell delay plus the mean of the arrivals.
        """
        if self.wireResistance is None: return
        position = instance.getAbutmentBox().getCenter()
        load     = self.getLoad( sinkPlugs, position )
        timing   = self.library.getCell( str(instance.getMasterCell().getName()) )
        delay    = ( timing.lookup('cell_rise', self.maxSlew, load)
                   + timing.lookup('cell_fall', self.maxSlew, load) ) / 2
        arrivals = self.getArrivals( sinkPlugs, position )
        self.latencies[ instance.getId() ] = delay + sum(arrivals) / len(arrivals)

    def __call__ ( self, sinkPlugs, position=None ):
        load = self.getLoad( sinkPlugs, position )
        for cell, timing in self.buffers:
            output = getBufferIos( cell )[1].getName()
            if load > timing.maxCap.get( str(output), load ): continue
            slew = max( timing.lookup('rise_transition', self.maxSlew, load)
                      , timing.lookup('fall_transition', self.maxSlew, load) )
            if slew <= self.maxSlew:
                return cell
        return self.buffers[-1][0]


def buildClockTree ( etesian, net, bufferCells, sliceHeight, maxSinks, maxRows=10 ):
    """
    Build, after placement, a buffer tree on ``net``. The tree is built
    bottom-up: the sinks of ``net`` are partitioned in groups (see
    partitionSinks()), each group gets its own buffer (see insertBuffer()),
    then the inputs of those buffers, which are now the sinks of ``net``,
    are partitioned in turn, until ``net`` itself has no more sinks than
    allowed at that level.

    ``bufferCells`` and ``maxSinks`` are either one value, or a list giving
    the value of each level, starting from the one closest to the original
    sinks (level 0), the last value being used for all the upper levels.
    ``bufferCells`` may also be a callable (for instance a BufferSizer),
    returning the buffer cell of a group of sinks.

    Buffers are named "cts_<net>_<level>_<index>".

    Return the list of the levels, each level being the list of its
    buffer instances.
    """
    def atLevel ( values, level ):
        if not isinstance(values, (list, tuple)): return values
        return values[ min(level, len(values)-1) ]

    for value in (maxSinks if isinstance(maxSinks, (list, tuple)) else [maxSinks]):
        if value < 2:
            raise ErrorMessage( 1, 'cts.buildClockTree(): maxSinks must be at least 2 (got {}).'.format( value ))
    levels = []
    while True:
        sinks      = getSinks( net )
        levelSinks = atLevel( maxSinks, len(levels) )
        if len(sinks) <= levelSinks: break
        level = []
        for group in partitionSinks( sinks, levelSinks ):
            name = 'cts_{}_{}_{}'.format( net.getName(), len(levels), len(level) )
            position = None
            if isinstance(bufferCells, BufferSizer):
                position   = bufferCells.getPosition( group )
                bufferCell = bufferCells( group, position )
            elif callable(bufferCells):
                bufferCell = bufferCells( group )
            else:
                bufferCell = atLevel( bufferCells, len(levels) )
            instance, outputNet = insertBuffer( etesian, net, group, bufferCell
                                              , name, sliceHeight, maxRows, position )
            if instance is None:
                raise ErrorMessage( 1, [ 'cts.buildClockTree(): No room left to place buffer "{}" of "{}".' \
                                         .format( name, net.getName() )
                                       , 'Searched {} rows around the sinks center, please increase free space.' \
                                         .format( maxRows ) ] )
            if isinstance(bufferCells, BufferSizer):
                bufferCells.setLatency( instance, group )
            level.append( instance )
        levels.append( level )
    return levels

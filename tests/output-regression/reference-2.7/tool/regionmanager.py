##################### generated by xml-casa (v2) from regionmanager.xml #############
##################### 75d5aaae93fa535a6d33fe8b20b02d65 ##############################
from __future__ import absolute_import 
from .__casac__ import regionmanager as _regionmanager
from .platform import str_encode as _str_ec
from .platform import str_decode as _str_dc
from .platform import dict_encode as _dict_ec
from .platform import dict_decode as _dict_dc
from .platform import encode as _any_ec
from .platform import decode as _any_dc
from .typecheck import validator as _pc
from .coercetype import coerce as _coerce


class regionmanager:
    ### self
    def __init__(self, *args, **kwargs):
        """This is the only regionmanager constructor.  It should generally be
        unnecessary for you to make one as there is little state in a
        regionmanager (you can set a Coordinate System with
        setcoordinates); the
        default regionmanager {stf rg} should be all you need.
        """
        self._swigobj = kwargs.get('swig_object',None)
        if self._swigobj is None:
            self._swigobj = _regionmanager()

    def absreltype(self, absrelvalue=int(0)):
        """This function is not intended for general user use.
        
        Regions may be specified with coordinates which are absolute or
        relative.  This function converts the integer code defining the
        absolute/relative type of the coordinates (which is stored in the
        region) into a string (maybe for printing purposes).
        
        The different types are
        
        
        Integer     String      Description
        1            abs        Absolute coordinate
        2            relref     Relative reference pixel
        3            relcen     Relative to center of image
        4            reldir     Relative to some direction
        
        
        """
        schema = {'absrelvalue': {'type': 'cInt'}}
        doc = {'absrelvalue': absrelvalue}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _absreltype_result = _str_dc(self._swigobj.absreltype(_pc.document['absrelvalue']))
        return _absreltype_result

    def box(self, blc=[ float(0) ], trc=[ float(-1) ], inc=[ float(1) ], absrel='abs', frac=False, comment=''):
        """This function creates a multi-dimensional pixel box region.  The box is
        specified by a bottom-left corner, and top-right corner and an increment
        (or stride).  Pixel coordinates are considered to run from 1 at the
        bottom left corner of the image to the image shape at the top-right
        corner of the image.
        
        You can specify whether the coordinates are given as pixel coordinates
        ({stfaf frac=F}) or fractions of the image shape ({stfaf frac=T}).
        Absolute fractions are in the range [0,1].
        
        You can also specify whether the coordinates are given as absolute
        coordinates ({stfaf absrel='abs'}) or relative to the reference pixel
        ({stfaf absrel='relref'}) or relative to the center of the image
        ({stfaf absrel='relcen'}).
        """
        schema = {'blc': {'type': 'cFloatVec', 'coerce': [_coerce.to_list,_coerce.to_floatvec]}, 'trc': {'type': 'cFloatVec', 'coerce': [_coerce.to_list,_coerce.to_floatvec]}, 'inc': {'type': 'cFloatVec', 'coerce': [_coerce.to_list,_coerce.to_floatvec]}, 'absrel': {'type': 'cStr'}, 'frac': {'type': 'cBool'}, 'comment': {'type': 'cStr'}}
        doc = {'blc': blc, 'trc': trc, 'inc': inc, 'absrel': absrel, 'frac': frac, 'comment': comment}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _box_result = _dict_dc(self._swigobj.box(_pc.document['blc'], _pc.document['trc'], _pc.document['inc'], _str_ec(_pc.document['absrel']), _pc.document['frac'], _str_ec(_pc.document['comment'])))
        return _box_result

    def frombcs(self, csys={ }, shape=[ int(0) ], box='', chans='', stokes='', stokescontrol='a', region=[ ]):
        """This function creates a multi-dimensional world coordinate region based
        on box, chans, stokes inputs familiar from image analysis tasks. It is
        being introduced as a temporary means of refactoring some python level
        task code into C++. However, if users find it to have value, its existence
        can be permanent.
        """
        schema = {'csys': {'type': 'cDict'}, 'shape': {'type': 'cIntVec', 'coerce': [_coerce.to_list,_coerce.to_intvec]}, 'box': {'type': 'cStr'}, 'chans': {'type': 'cStr'}, 'stokes': {'type': 'cStr'}, 'stokescontrol': {'type': 'cStr'}, 'region': {'type': 'cVariant'}}
        doc = {'csys': csys, 'shape': shape, 'box': box, 'chans': chans, 'stokes': stokes, 'stokescontrol': stokescontrol, 'region': region}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _frombcs_result = _dict_dc(self._swigobj.frombcs(_dict_ec(_pc.document['csys']), _pc.document['shape'], _str_ec(_pc.document['box']), _str_ec(_pc.document['chans']), _str_ec(_pc.document['stokes']), _str_ec(_pc.document['stokescontrol']), _any_ec(_pc.document['region'])))
        return _frombcs_result

    def complement(self, region=[ ], comment=''):
        """This function (short-hand name {tt comp}) creates the complement of
        a world region(s).
        
        The region parameter can be a single region record defining a simple
        or complex region or it can contain several region records in a
        Python dictionary.  If multiple regions are given then the union of
        this set of regions is taken first, and the complement is found from
        the union.
        
        NOTE: ia.statistics() is UNABLE to handle complement regions in CASA yet.
        """
        schema = {'region': {'type': 'cVariant'}, 'comment': {'type': 'cStr'}}
        doc = {'region': region, 'comment': comment}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _complement_result = _dict_dc(self._swigobj.complement(_any_ec(_pc.document['region']), _str_ec(_pc.document['comment'])))
        return _complement_result

    def concatenation(self, box=[ ], regions=[ ], comment=''):
        """This function (short-hand name {tt concat}) creates a region which is
        the concatenation along a new axis of the given world regions.
        
        This function is similar to the
        extension function.  The
        {stfaf concatenation} function allows you to take many world regions,
        and concatenate them along one axis (rather than take one region and
        extend it along many axes which is what function {stff extension}
        does).
        
        For example, you may have generated a different polygonal region for
        each spectral pixel of a spectral-line cube and you wish to concatenate them
        together to form the overall region for use in a deconvolution
        application.
        
        The axis to concatenate along is specified as a 1-dimensional world box.
        The shape of the 1D box must contain as many pixels (although you
        don't have to specify it in pixels) as there are regions
        to concatenate.
        
        Because this function is most likely to be used in a script, the
        interface takes a record containing {stff region} records, Python
        dictionaries, as there might be a lot of them.
        """
        schema = {'box': {'type': 'cVariant'}, 'regions': {'type': 'cVariant'}, 'comment': {'type': 'cStr'}}
        doc = {'box': box, 'regions': regions, 'comment': comment}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _concatenation_result = _dict_dc(self._swigobj.concatenation(_any_ec(_pc.document['box']), _any_ec(_pc.document['regions']), _str_ec(_pc.document['comment'])))
        return _concatenation_result

    def deletefromtable(self, tablename='', regionname=''):
        """This function deletes a region stored in an casa  Table.
        
        For the {stfaf tablename} argument,
        
        you have to give  the name of an existing
        CASA table on disk (any kind of table).
        
        You specify the name of the region with the {stfaf regionname}
        arguments.  If you set {stfaf regionname=''} then nothing is done.  The names of all the regions stored in a Table can be found
        with the function
        namesintable.
        """
        schema = {'tablename': {'type': 'cStr'}, 'regionname': {'type': 'cStr'}}
        doc = {'tablename': tablename, 'regionname': regionname}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _deletefromtable_result = self._swigobj.deletefromtable(_str_ec(_pc.document['tablename']), _str_ec(_pc.document['regionname']))
        return _deletefromtable_result

    def difference(self, region1={ }, region2={ }, comment=''):
        """This function (short-hand name {stff diff}) creates
        a region which is the difference of two world regions.  The order
        of the regions is important.
        
        The difference consists of all pixels masked-on in the first
        region and not masked-on in the second region.
        """
        schema = {'region1': {'type': 'cDict'}, 'region2': {'type': 'cDict'}, 'comment': {'type': 'cStr'}}
        doc = {'region1': region1, 'region2': region2, 'comment': comment}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _difference_result = _dict_dc(self._swigobj.difference(_dict_ec(_pc.document['region1']), _dict_ec(_pc.document['region2']), _str_ec(_pc.document['comment'])))
        return _difference_result

    def done(self):
        """This function destroys the contents of the {stf regionmanager} tool
        (including its GUI).  The tool still exists as a Glish variable, but
        it is no longer a Regionmanager ! You are unlikely to need this
        function.
        
        """
        _done_result = self._swigobj.done()
        return _done_result

    def selectedchannels(self, specification='', shape=[ int(0) ]):
        """This method returns all the selected zero-based channel numbers from the specified string within the image.
        
        
        
        
        """
        schema = {'specification': {'type': 'cStr'}, 'shape': {'type': 'cIntVec', 'coerce': [_coerce.to_list,_coerce.to_intvec]}}
        doc = {'specification': specification, 'shape': shape}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _selectedchannels_result = self._swigobj.selectedchannels(_str_ec(_pc.document['specification']), _pc.document['shape'])
        return _selectedchannels_result

    def fromtextfile(self, filename='', shape=[ int(0) ], csys={ }):
        """This function reads a text file containing region descriptions and
        converts it to a python dictionary.
        
        
        """
        schema = {'filename': {'type': 'cReqPath', 'coerce': _coerce.expand_path}, 'shape': {'type': 'cIntVec', 'coerce': [_coerce.to_list,_coerce.to_intvec]}, 'csys': {'type': 'cDict'}}
        doc = {'filename': filename, 'shape': shape, 'csys': csys}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _fromtextfile_result = _dict_dc(self._swigobj.fromtextfile(_str_ec(_pc.document['filename']), _pc.document['shape'], _dict_ec(_pc.document['csys'])))
        return _fromtextfile_result

    def fromtext(self, text='', shape=[ int(1) ], csys={ }):
        """This function reads a region region text descriptions and
        converts it to a python region dictionary.
        
        
        """
        schema = {'text': {'type': 'cStr'}, 'shape': {'type': 'cIntVec', 'coerce': [_coerce.to_list,_coerce.to_intvec]}, 'csys': {'type': 'cDict'}}
        doc = {'text': text, 'shape': shape, 'csys': csys}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _fromtext_result = _dict_dc(self._swigobj.fromtext(_str_ec(_pc.document['text']), _pc.document['shape'], _dict_ec(_pc.document['csys'])))
        return _fromtext_result

    def fromfiletorecord(self, filename='', verbose=True, regionname=''):
        """This function reads files containing ImageRegion objects and turns them
        into Region Records.
        
        The intended use for this method is to read the file saved by the casa
        viewer and turn the files contents into regions that are usabla by the
        image analysis tool.
        """
        schema = {'filename': {'type': 'cStr'}, 'verbose': {'type': 'cBool'}, 'regionname': {'type': 'cStr'}}
        doc = {'filename': filename, 'verbose': verbose, 'regionname': regionname}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _fromfiletorecord_result = _dict_dc(self._swigobj.fromfiletorecord(_str_ec(_pc.document['filename']), _pc.document['verbose'], _str_ec(_pc.document['regionname'])))
        return _fromfiletorecord_result

    def tofile(self, filename='', region={ }):
        """This function is to store a region created by the regionmanager in a disk file for future use
        """
        schema = {'filename': {'type': 'cStr'}, 'region': {'type': 'cDict'}}
        doc = {'filename': filename, 'region': region}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _tofile_result = self._swigobj.tofile(_str_ec(_pc.document['filename']), _dict_ec(_pc.document['region']))
        return _tofile_result

    def fromrecordtotable(self, tablename='', regionname=[ ], regionrec={ }, asmask=False, verbose=True):
        """This function saves regions into an casa Table
        For the {stfaf tablename} argument the user should be the name of an existing
        aipspp Table on disk (any kind of table).
        
        If the parameter {tt asmask} is {tt True} then the table has to be an image table.
        A mask makes sense with an image only.
        
        
        You can specify the name the region will have ({stfaf
        regionname}) when it is saved in the Table.  If you don't specify this,
        a digit based name is assigned to it or if specify a name that already
        exists a new one will be generated which is close but different. The
        function returns you the name the region is assigned
        """
        schema = {'tablename': {'type': 'cStr'}, 'regionname': {'type': 'cVariant'}, 'regionrec': {'type': 'cDict'}, 'asmask': {'type': 'cBool'}, 'verbose': {'type': 'cBool'}}
        doc = {'tablename': tablename, 'regionname': regionname, 'regionrec': regionrec, 'asmask': asmask, 'verbose': verbose}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _fromrecordtotable_result = _str_dc(self._swigobj.fromrecordtotable(_str_ec(_pc.document['tablename']), _any_ec(_pc.document['regionname']), _dict_ec(_pc.document['regionrec']), _pc.document['asmask'], _pc.document['verbose']))
        return _fromrecordtotable_result

    def fromtabletorecord(self, tablename='', regionname=[ ], verbose=True):
        """This function restores a region from an aipspp Table
        to the global name space.
        
        For the {stfaf tablename} argument, you can specify an
        image tool, a table tool,
        or a string.  If you give a string, it should be the name of an existing
        aipspp table on disk (any kind of table).
        
        If {stfaf numberfields} is F, then the field names of the
        record are the same as they are in the Table.  Otherwise,
        the regions are put into numbered fields (the field
        names could be anything).
        
        You can use the function
        namesintable to find out the
        names of the regions in the Table.
        """
        schema = {'tablename': {'type': 'cStr'}, 'regionname': {'type': 'cVariant'}, 'verbose': {'type': 'cBool'}}
        doc = {'tablename': tablename, 'regionname': regionname, 'verbose': verbose}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _fromtabletorecord_result = _dict_dc(self._swigobj.fromtabletorecord(_str_ec(_pc.document['tablename']), _any_ec(_pc.document['regionname']), _pc.document['verbose']))
        return _fromtabletorecord_result

    def intersection(self, regions=[ ], comment=''):
        """This function (short-hand name {stff int}) creates a region which is
        the intersection of the given world regions.   The input regions can
        themselves be compound regions (such as the union or intersection etc).
        The input regions must be provided as a Python dictionary of regions
        (see examples).
        """
        schema = {'regions': {'type': 'cVariant'}, 'comment': {'type': 'cStr'}}
        doc = {'regions': regions, 'comment': comment}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _intersection_result = _dict_dc(self._swigobj.intersection(_any_ec(_pc.document['regions']), _str_ec(_pc.document['comment'])))
        return _intersection_result

    def ispixelregion(self, region={ }):
        """NOT IMPLEMENTED IN CASA
        
        This function returns T if the region is a pixel region.
        For any other glish variable it returns F.
        """
        schema = {'region': {'type': 'cDict'}}
        doc = {'region': region}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _ispixelregion_result = self._swigobj.ispixelregion(_dict_ec(_pc.document['region']))
        return _ispixelregion_result

    def isworldregion(self, region={ }):
        """NOT IMPLEMENTED IN CASA
        
        This function returns T if the region is a world region.
        For any other glish variable it returns F.
        """
        schema = {'region': {'type': 'cDict'}}
        doc = {'region': region}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _isworldregion_result = self._swigobj.isworldregion(_dict_ec(_pc.document['region']))
        return _isworldregion_result

    def namesintable(self, tablename=''):
        """This function returns the names of regions stored in an CASA Table.
        
        For the {stfaf tablename} argument, you can specify a string; it should be the name of an existing
        aipspp table on disk (any kind of table).
        """
        schema = {'tablename': {'type': 'cStr'}}
        doc = {'tablename': tablename}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _namesintable_result = [_str_dc(_x) for _x in self._swigobj.namesintable(_str_ec(_pc.document['tablename']))]
        return _namesintable_result

    def setcoordinates(self, csys={ }):
        """This function allows you to (re)set the default Coordinate System
        used by the functions that make world regions.  If you don't specifiy a
        Coordinate System when you make the world region, the default Coordinate
        System, if there is one, is used.   The Coordinate System is
        stored in a {stf coordinates} tool and is created with
        the coordsys toolfunction.
        
        Normally, the world region creating functions like
        wbox and
        wpolygon will issue a message
        each time the private Coordinate System is used.  However, if you set
        {stfaf verbose=F} then this will not occur.
        
        
        """
        schema = {'csys': {'type': 'cDict'}}
        doc = {'csys': csys}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _setcoordinates_result = self._swigobj.setcoordinates(_dict_ec(_pc.document['csys']))
        return _setcoordinates_result

    def makeunion(self, regions=[ ], comment=''):
        """This function takes a minimum of two world regions and creates a region which
        is the union of the given regions.  The input regions can themselves be
        compound regions (such as the union or intersection etc).   The input
        regions must be a Pythion dictionary of at leat two regions
        (see examples).
        """
        schema = {'regions': {'type': 'cVariant'}, 'comment': {'type': 'cStr'}}
        doc = {'regions': regions, 'comment': comment}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _makeunion_result = _dict_dc(self._swigobj.makeunion(_any_ec(_pc.document['regions']), _str_ec(_pc.document['comment'])))
        return _makeunion_result

    def wbox(self, blc=[ ], trc=[ ], pixelaxes=[ int(-1) ], csys={ }, absrel='abs', comment=''):
        """This function creates a multi-dimensional world box region; the
        corners of the box are specified in world coordinates.  However, the box
        is not a true world volume in that its sides do not follow world
        contours.  Its sides are parallel to the pixel axes.  If you are in a
        region of high world coordinate contour non-linearity (e.g.  near a
        pole), you are probably better off using a world polygon.
        
        The box is specified by a bottom-left corner, and a top-right corner.
        The coordinates are given as quantities, and you can give a vector of
        quantities (e.g.  {cf blc = qa.quantity("1rad 20deg")} or a
        quantity of a vector (e.g.{cf blc = qa.quantity([10,30], 'rad')}).
        
        You can specify whether the coordinates are given as absolute coordinates
        ({stfaf absrel='abs'}) or relative to the reference pixel ({stfaf
        absrel='relref'}) or relative to the center of the image ({stfaf
        absrel='relcen'}).  You can specify this for each axis (the same for the
        blc and trc).   If you specify less values than the number of
        values in {stfaf blc} or {stfaf trc} then the last value you
        did specify is used as the default for all higher numbered axes
        (e.g. {stfaf absrel='relref'} means {stfaf absrel="relref relref"}
        for two axes).
        
        You specify which pixel axes in the image the {stfaf blc} and {stfaf
        trc} vector refer to with the {stfaf pixelaxes} argument.  If you
        don't, it defaults to [0,1,2,...].  This specification is an important
        part of world regions.
        
        You must also specify the Coordinate System with the {stfaf csys}
        argument.  The Coordinate System is encapsulated in a {stfaf coordinates}
        tool and can be recovered from an image with the
        coordsys tool function.  You can
        also set a default Coordinate System in the regionmanager with the
        setcoordinates
        function.
        
        In the regionmanager we have defined units `pix' and `frac'; these are
        then known to the quanta system.  This means
        that you can effectively define a pixel box (except for the stride
        capability) as a world box with most of the advantages of world regions
        (can be used for compound regions).  However, it is still not very
        portable to other images because the coordinates are pixel based,
        not world based.
        
        Note that the need to deal with the {stfaf pixelaxes} and {stfaf csys}
        is hidden from you when using the gui
        interface of the regionmanager.
        """
        schema = {'blc': {'type': 'cVariant'}, 'trc': {'type': 'cVariant'}, 'pixelaxes': {'type': 'cIntVec', 'coerce': [_coerce.to_list,_coerce.to_intvec]}, 'csys': {'type': 'cDict'}, 'absrel': {'type': 'cStr'}, 'comment': {'type': 'cStr'}}
        doc = {'blc': blc, 'trc': trc, 'pixelaxes': pixelaxes, 'csys': csys, 'absrel': absrel, 'comment': comment}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _wbox_result = _dict_dc(self._swigobj.wbox(_any_ec(_pc.document['blc']), _any_ec(_pc.document['trc']), _pc.document['pixelaxes'], _dict_ec(_pc.document['csys']), _str_ec(_pc.document['absrel']), _str_ec(_pc.document['comment'])))
        return _wbox_result

    def wpolygon(self, x=[ ], y=[ ], pixelaxes=[ int(-1) ], csys={ }, absrel='abs', comment=''):
        """This function (short-hand name {stff wpoly}) creates a 2D world
        polygon region.  The polygon is specified by an {stfaf x} and a {stfaf y}
        vector.  These must be quantities of a vector (the
        world box function
        allows both
        quantities of vectors and vectors of quantities).  This means that the
        units are common to all elements of each vector.  Thus, {cf
        qa.quantity([1,2,3],'rad')} (a quantity of a vector) is different from
        {cf qa.quantity("1rad 2rad 3rad")} (a vector of quantities) although
        the information that they specify is the same.
        
        You specify which pixel axes in the image the {stfaf x} and {stfaf
        y} vectors pertain to with the {stfaf pixelaxes} argument.  If you don't,
        it defaults to [0,1].  This specification is an important part of
        world regions.
        
        You can specify whether the {stfaf x} and {stfaf y} vector coordinates are
        given as absolute coordinates ({stfaf absrel='abs'}) or relative to the
        reference pixel ({stfaf absrel='relref'}) or relative to the center of the
        image ({stfaf absrel='relcen'}).  This argument applies to both the axes
        of the polygon.
        
        You must also specify the Coordinate System with the {stfaf csys}
        argument.  The Coordinate System is encapsulated in a {stfaf coordinates}
        tool and can be recovered from an image with the
        coordsys function.  You can
        also set a default Coordinate System in the Regionmanager with the
        setcoordinates
        function.
        
        In the regionmanager we have defined units `pix' and `frac'; these are
        then known to the quanta system.  This means
        that you can effectively define a pixel box (except for the stride
        capability) as a world box with most of the advantages of world regions
        (can be used for compound regions).  However, it is still not very
        portable to other images because the coordinates are pixel based,
        not world based.
        
        Note that the need to deal with the {stfaf pixelaxes} and {stfaf csys}
        is hidden from you when using the gui
        interface of the regionmanager.
        """
        schema = {'x': {'type': 'cVariant'}, 'y': {'type': 'cVariant'}, 'pixelaxes': {'type': 'cIntVec', 'coerce': [_coerce.to_list,_coerce.to_intvec]}, 'csys': {'type': 'cDict'}, 'absrel': {'type': 'cStr'}, 'comment': {'type': 'cStr'}}
        doc = {'x': x, 'y': y, 'pixelaxes': pixelaxes, 'csys': csys, 'absrel': absrel, 'comment': comment}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _wpolygon_result = _dict_dc(self._swigobj.wpolygon(_any_ec(_pc.document['x']), _any_ec(_pc.document['y']), _pc.document['pixelaxes'], _dict_ec(_pc.document['csys']), _str_ec(_pc.document['absrel']), _str_ec(_pc.document['comment'])))
        return _wpolygon_result


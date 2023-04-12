classdef Device < types.core.NWBContainer & types.untyped.GroupClass
% DEVICE Metadata about a data acquisition device, e.g., recording system, electrode, microscope.


% PROPERTIES
properties
    description; % Description of the device (e.g., model, firmware version, processing software version, etc.) as free-form text.
    manufacturer; % The name of the manufacturer of the device.
end

methods
    function obj = Device(varargin)
        % DEVICE Constructor for Device
        %     obj = DEVICE(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % description = char
        % manufacturer = char
        obj = obj@types.core.NWBContainer(varargin{:});
        
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'description',[]);
        addParameter(p, 'manufacturer',[]);
        misc.parseSkipInvalidName(p, varargin);
        obj.description = p.Results.description;
        obj.manufacturer = p.Results.manufacturer;
        if strcmp(class(obj), 'types.core.Device')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    function obj = set.description(obj, val)
        obj.description = obj.validate_description(val);
    end
    function obj = set.manufacturer(obj, val)
        obj.manufacturer = obj.validate_manufacturer(val);
    end
    %% VALIDATORS
    
    function val = validate_description(obj, val)
        val = types.util.checkDtype('description', 'char', val);
    end
    function val = validate_manufacturer(obj, val)
        val = types.util.checkDtype('manufacturer', 'char', val);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.NWBContainer(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        if ~isempty(obj.description)
            io.writeAttribute(fid, [fullpath '/description'], obj.description);
        end
        if ~isempty(obj.manufacturer)
            io.writeAttribute(fid, [fullpath '/manufacturer'], obj.manufacturer);
        end
    end
end

end
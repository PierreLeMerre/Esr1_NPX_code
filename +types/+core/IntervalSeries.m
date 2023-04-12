classdef IntervalSeries < types.core.TimeSeries & types.untyped.GroupClass
% INTERVALSERIES Stores intervals of data. The timestamps field stores the beginning and end of intervals. The data field stores whether the interval just started (>0 value) or ended (<0 value). Different interval types can be represented in the same series by using multiple key values (eg, 1 for feature A, 2 for feature B, 3 for feature C, etc). The field data stores an 8-bit integer. This is largely an alias of a standard TimeSeries but that is identifiable as representing time intervals in a machine-readable way.



methods
    function obj = IntervalSeries(varargin)
        % INTERVALSERIES Constructor for IntervalSeries
        %     obj = INTERVALSERIES(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        varargin = [{'data_resolution' types.util.correctType(-1.0, 'float32') 'data_unit' 'n/a'} varargin];
        obj = obj@types.core.TimeSeries(varargin{:});
        if strcmp(class(obj), 'types.core.IntervalSeries')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    
    %% VALIDATORS
    
    function val = validate_data(obj, val)
        val = types.util.checkDtype('data', 'int8', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = val.dims;
        else
            valsz = size(val);
        end
        validshapes = {[Inf]};
        types.util.checkDims(valsz, validshapes);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.TimeSeries(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
    end
end

end
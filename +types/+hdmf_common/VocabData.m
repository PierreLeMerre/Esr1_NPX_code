classdef VocabData < types.hdmf_common.VectorData & types.untyped.DatasetClass
% VOCABDATA Data that come from a controlled vocabulary of text values. A data value of i corresponds to the i-th element in the 'vocabulary' array attribute.


% PROPERTIES
properties
    vocabulary; % The available items in the controlled vocabulary.
end

methods
    function obj = VocabData(varargin)
        % VOCABDATA Constructor for VocabData
        %     obj = VOCABDATA(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % vocabulary = char
        obj = obj@types.hdmf_common.VectorData(varargin{:});
        
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'vocabulary',[]);
        misc.parseSkipInvalidName(p, varargin);
        obj.vocabulary = p.Results.vocabulary;
        if strcmp(class(obj), 'types.hdmf_common.VocabData')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    function obj = set.vocabulary(obj, val)
        obj.vocabulary = obj.validate_vocabulary(val);
    end
    %% VALIDATORS
    
    function val = validate_data(obj, val)
        val = types.util.checkDtype('data', 'uint8', val);
    end
    function val = validate_description(obj, val)
        val = types.util.checkDtype('description', 'char', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = val.dims;
        elseif istable(val)
            valsz = height(val);
        elseif ischar(val)
            valsz = size(val, 1);
        else
            valsz = size(val);
        end
        validshapes = {[1]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_vocabulary(obj, val)
        val = types.util.checkDtype('vocabulary', 'char', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = val.dims;
        elseif istable(val)
            valsz = height(val);
        elseif ischar(val)
            valsz = size(val, 1);
        else
            valsz = size(val);
        end
        validshapes = {[1]};
        types.util.checkDims(valsz, validshapes);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.hdmf_common.VectorData(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        if ~isempty(obj.vocabulary)
            io.writeAttribute(fid, [fullpath '/vocabulary'], obj.vocabulary);
        else
            error('Property `vocabulary` is required in `%s`.', fullpath);
        end
    end
end

end
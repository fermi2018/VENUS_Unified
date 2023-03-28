% clear
% load
s = whos;
for i = 1:length(s)
    if strcmp(s(i).class,'mp')
        name = s(i).name;
        assignin('base', name, double(eval('base', name)));
    elseif (strcmp(s(i).class,'struct')==1 && max(s(i).size)==1)
        name=s(i).name;
        fn = fieldnames(eval(name));
        for k=1:numel(fn)
            if( isnumeric(getfield(eval(name),fn{k})) )
                % 
                assignin('base',name,setfield(eval(name),fn{k},double(eval([name,'.',fn{k}]))));
            elseif ( iscell(getfield(eval(name),fn{k})) && isnumeric(eval([name,'.',fn{k},'{1}'])))
%                 keyboard
                for ic = 1:length(eval([name,'.',fn{k}]))
                    nameIc=[name,'.',fn{k},'{',num2str(ic),'}'];
                    v=double(eval(nameIc));
                    evalin('base', [nameIc,'=v']);
                end
            end
        end
    elseif strcmp(s(i).class,'struct')==1
        name=s(i).name;
        fnIs = fieldnames(eval(name));
        for is = 1:length(eval(name))
            nameIs=[name,'(',num2str(is),')'];
            for kIs=1:numel(fnIs)
                if( isnumeric(getfield(eval(nameIs),fnIs{kIs})) )
                    % do stuff
%                     keyboard
                    
%                     assignin('base',nameIs,setfield(eval(nameIs),fnIs{kIs},double(eval([nameIs,'.',fnIs{kIs}]))));
                    geom(is).(fnIs{kIs})=double(geom(is).(fnIs{kIs}));
%                      a=deal(setfield(eval(nameIs),fnIs{kIs},double(eval([nameIs,'.',fnIs{kIs}]))));
%                         assignin('base',nameIs,a)
                end
                
            end
        end
    end
end

multiprecision=0;
clear is fnIs kIs name nameIs i s name


function [n_dba,n_props,q0,acc,q_prime]=read_spec(filename)
n_dba= 0;
n_props= 0;
q0= 0;
acc= 0;

fid = fopen(filename);
tline = fgetl(fid);
while ischar(tline)
    str= split(tline,["=","#",","]);
    if(strcmp(str{1}, 'Name'))
        disp(tline)
    end
    if(strcmp(str{1}, 'AP'))
        disp(tline)
    end
    if(strcmp(str{1}, 'nAP'))
        n_props = 2^(str2double(str{2}));
    end
    if(strcmp(str{1}, 'nNodes'))
        n_dba= str2double(str{2});
    end
    if(strcmp(str{1}, 'init'))
        q0= str2double(str{2});
    end
    if(strcmp(str{1}, 'acc'))
        acc= str2double(str{2});
    end
    tline= fgetl(fid);
end
q_prime= zeros(n_dba, n_props);

frewind(fid);
tline = fgetl(fid);
while ischar(tline)
    str= split(tline,["=","#",","]);
    if(numel(str) == 3)
        q_prime(uint8(str2num(str{1}))+1, ...
            uint8(str2num(str{2}))+1)= uint8(str2num(str{3}))+1;
    end
    tline= fgetl(fid);
end

fclose(fid);

end
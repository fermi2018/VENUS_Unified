function n=extract_n(Nline,chs)
%         chs='diam=';
         ise=length(chs)+1;
         ipu=findstr(Nline,chs);
         [chr,remainder]=strtok(Nline(ipu:end));
         n=str2num(chr(ise:end));

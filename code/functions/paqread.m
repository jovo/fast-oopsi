function [data]=paqread(filename, varargin);
% function [data]=paqread(filename, varargin);
%
% Read binary files acquired with PackIO for use in DaqViewer
% Return recording of each channel as a vector array
% Based on original program from JJ
% 1st input = name of file to open and read
% 2nd input = read mode: either 'channels' or 'info'
% Optional 3rd input = channels, vector list of channel numbers to read

% Adam Packer
% January 17th, 2007


fid=fopen(filename);

rate=fread(fid,1,'float32','b');
numchans=fread(fid,1,'float32','b');
for i=1:numchans;
    number_of_characters=fread(fid,1,'float32','b');
    channelname{i}=[];
    for j=1:number_of_characters
        channelname{i}=[channelname{i}, strrep(fread(fid,1,'float32=>char', 'b'),' ','')];
    end
end

[pathstr, name, ext, versn] = fileparts(filename);
if strcmp(ext,'.paq')
    for k=1:numchans
        number_of_characters=fread(fid,1,'float32','b');
        HWchan{k}=[];
        for m=1:number_of_characters
            HWchan{k}=[HWchan{k}, strrep(strrep(fread(fid,1,'float32=>char','b'),' ',''),'/','_')];
        end
    end
    
    for n=1:numchans
        number_of_characters=fread(fid,1,'float32','b');
        units{n}=[];
        for q=1:number_of_characters
            units{n}=[units{n}, strrep(fread(fid,1,'float32=>char','b'),' ','')];
        end
    end
end

fposition=ftell(fid);
[pathstr, name, ext, versn] = fileparts(filename);
dirinfo=dir(pathstr);
for n=1:size(dirinfo);
    if strcmp(dirinfo(n).name,[name ext]);
        filesize=dirinfo(n).bytes;
    end
end

switch lower(varargin{1})
    case 'channels' % DATA = daqread(FILENAME);
        includedchannels=varargin{2};
        try
            [rangequest,qmquestion]=ephysdefaults;
        catch
            rangequest=0;
            qmquestion=0;
        end

        if rangequest==1;
            prompt1=['Filesize(MB) = ',...
                sprintf('%-.2f\n',filesize/10^6),...
                '1sec data (MB) = ',...
                sprintf('%-.2f\n',4*rate*numchans/10^6),...
                'Filesize(sec) = ',...
                sprintf('%-.2f\n',(filesize/10^6)/(4*rate*numchans/10^6)),...
                'Read from time(sec)'];
            prompt={prompt1 'Stop reading at time(sec) [Enter zero for all]'};
            dlg_title=['Select range to import:'];
            num_lines=1;
            def={'0','0'};
            answer=inputdlg(prompt,dlg_title,num_lines,def);
            starttime=str2num(cell2mat(answer(1)));
            stoptime=str2num(cell2mat(answer(2)));
            startsample=starttime*rate;
            stopsample=stoptime*rate;
            startbyte=startsample*numchans;
            stopbyte=stopsample*numchans;
            if (stopbyte > filesize) ||(stopbyte==0)
                stopbyte=filesize;
            end

            if qmquestion==0;
                button='Quick';
            elseif qmquestion==1;
                button = questdlg('Read quickly OR Open more samples?','Reading style','Quick','More','Quick');
            end

            if strcmp(button,'More')        % READ BY FOR LOOP THROUGH DATA
                waithandle = waitbar(0,'Reading from PackIO file');%for user
                for i=1:length(includedchannels)
                    fseek(fid,fposition+startbyte*4+includedchannels(i)*4-4,'bof');
                    data(:,i)=fread(fid,(stopbyte-startbyte)/numchans,'*float32',4*(numchans-1),'b');
                    waithandle = waitbar(i/length(includedchannels),waithandle);
                end
                close(waithandle);
                fclose(fid);
            elseif strcmp(button,'Quick')   % READ BY AUTO COLUMN SORT FREAD
                fseek(fid,fposition+startbyte*4,'bof');
                data=fread(fid,[numchans,(stopbyte-startbyte)/numchans],'*float32','b');
                notchans=1:numchans;
                notchans(includedchannels)=[];
                data(notchans,:)=[];
                data=data';
                fclose(fid);
            end
        elseif rangequest==0;
            fseek(fid,fposition,'bof');
            data=fread(fid,[numchans,filesize/numchans],'*float32','b');
            notchans=1:numchans;
            notchans(includedchannels)=[];
            data(notchans,:)=[];
            data=data';
            fclose(fid);
        end

    case 'info' % DAQINFO = paqread(FILENAME, 'info');
        % Error if an invalid second argument is passed.
        if ~strcmp(lower(varargin{1}), 'info')
            fclose(fid);
            error('daq:daqread:invalidarg', 'Invalid second argument specified.');
        end
        % Read in the Object information and Event information.
        info.ObjInfo.SampleRate=rate;
        info.HwInfo='Probably acquired with PackIO';
        info.ObjInfo.SamplesAcquired=((filesize-fposition)/4)/numchans;
        for a=1:numchans;
            if strcmp(ext,'.paq')
                info.ObjInfo.Channel(a).HwChannel=HWchan{a};
                info.ObjInfo.Channel(a).Units=units{a};
            else
                info.ObjInfo.Channel(a).HwChannel=0;
                info.ObjInfo.Channel(a).Units='SeePackIO';
            end
            info.ObjInfo.Channel(a).ChannelName=channelname{a};
        end
        data = info;
        fclose(fid);
end

% m=daqmem;
% memory=m.AvailPhys;
% uncomment below to test low mem functionality even if you have high mem
% if memory > filesize;
% if memory < filesize;
%     prompt={'Read from byte #:' 'Stop reading at byte #:'};
%     dlg_title='Memory is limited.';
%     num_lines=1;
%     def={'',''};
%     answer=inputdlg(prompt,dlg_title,num_lines,def);
%     start=str2num(cell2mat(answer(1)));
%     stop=str2num(cell2mat(answer(2)));
%     fseek(fid,start,'cof');
%     temp=ftell(fid);
%     data=fread(fid,stop-start,'*float32', 'b');
%     fclose(fid);
%     data=reshape(data, numchans, size(data, 1)/numchans);
%     for i=1:numchans
%         assignin('base', channelname{i}, data(i, :));
%     end
% else
%     for i=1:numchans
%         assignin('base', channelname{i}, data(i, :));
%     end
% end
% Time=0:size(data, 2)-1;
% Time=Time/rate;
% assignin('base', 'Time', Time);

% READ BY MEMORY MAP
% close if making memmap
% poss offset for memmap'Offset',fposition*4,...
%poss format for memmap             [(stopbyte-startbyte)/numchans numchans] 'x'
%                             fclose(fid)
%             m = memmapfile(filename,...
%                 'Format','single',...
%                 'Writable',false);
%             data=m.Data(startbyte:stopbyte,includedchanne
%             ls);
function [] = writeGoICPCfgFile( msethresh,  distTransSize)

f_id = fopen('GoICP_config.cfg', 'w'); 

fprintf(f_id, ...
['# Config file for GO-ICP \n MSEThresh=%f \n rotMinX=-3.1416 \n rotMinY=-3.1416 \n rotMinZ=-3.1416 \n ',...
'rotWidth=6.2832 \n transMinX=-0.5 \n transMinY=-0.5 \n transMinZ=-0.5 \n transWidth=1.0 \n ',...
'trimFraction=0.1 \n distTransSize=%d \n distTransExpandFactor=2.0'], msethresh, distTransSize);

fclose(f_id);

end


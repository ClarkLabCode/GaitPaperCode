function [ newData ] = PrepareTurningData(newData, vfMin)
%Utility function to run all data preprocessing

[ newData ] = computeAnalyticSignal( newData );
[ newData ] = filterData( newData, [], [], 15, 3 );
[ newData ] = addSmoothedCameraFrameUpDown( newData );
[ newData ] = computeSwingStanceAmplitudesAndDurations( newData );
[ newData.yawExtremum ]= LabelYawExtrema(newData, 'findpeaks', vfMin);

end


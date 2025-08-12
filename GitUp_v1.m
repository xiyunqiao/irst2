clear all 
img_path = 'D:\Phd\03_Research_and_subject\00-Research\Infrared\data21';
files = dir(img_path);
frame_diff_path = 'D:\Phd\03_Research_and_subject\00-Research\写论文\03-GRSL_Temporal_Feature\My_method\data21_res_frame_diff';
output_path = 'D:\Phd\03_Research_and_subject\00-Research\写论文\03-GRSL_Temporal_Feature\My_method\data21_res_new_0729';
if ~exist(output_path, 'dir')
    % 如果文件夹不存在，则创建它
    mkdir(output_path);
end
if ~exist(frame_diff_path, 'dir')
    % 如果文件夹不存在，则创建它
    mkdir(frame_diff_path);
end
img_path_list = dir(fullfile(img_path,'*.bmp'));%获取该文件夹中所有bmp格式的图像
num = length(img_path_list);


files_name =sort_nat({img_path_list.name});%重新排序
for i = 1:num
    newname=fullfile(img_path, files_name{i});
    img1 = imread(newname);
    nn = ndims(img1);
    if nn==3
      img1= rgb2gray(img1);
    end
end


for i = 1:num
    t_thr = 5;
    if i>t_thr
        I_prev  = rgb2gray(imread(fullfile(img_path, files_name{i-t_thr})));
        I_mid   = rgb2gray(imread(fullfile(img_path, files_name{i})));
        I_next  = rgb2gray(imread(fullfile(img_path, files_name{i+t_thr})));
        [regis_tmp1, x1, y1] = registration_1(I_mid(:,:,1), I_prev(:,:,1));
        [regis_tmp2, x2, y2]  = registration_1(I_mid(:,:,1), I_next(:,:,1));
        x_max = ceil(max(x1, x2));
        y_max = ceil(max(y1, y2));
        D1 = (double(I_mid) - double(regis_tmp1));
        D2 = (double(regis_tmp2) - double(I_mid));

        D1(:, 1:x_max+1) = 0;
        D1(:, 256-x_max-1:256) = 0;
        D1(1:y_max+1, :) = 0;
        D1(256-y_max-1:256,:) = 0;

        D2(:, 1:x_max+1) = 0;
        D2(:, 256-x_max-1:256) = 0;
        D2(1:y_max+1, :) = 0;
        D2(256-y_max-1:256,:) = 0;

        frame_diff = mat2gray(max(D1, D2));


%         mask_mid = true(size(I_mid));  % 中间帧本身是完整的
%         common_mask = mask_prev & mask_next & mask_mid;
%         Diff_masked = frame_diff;
%         Diff_masked(~common_mask) = 0;
        file_str = fullfile(frame_diff_path, files_name{i});
        imwrite(mat2gray(frame_diff), file_str);
        histogram_res = candidate_mask_kourenke(frame_diff);
        [SCR_map, TEMPORAL_map] = Morphology_analysis(files_name{i}, histogram_res, img_path, t_thr, frame_diff);
        result = SCR_map.*TEMPORAL_map;
 
    else
        pic2 = imread(fullfile(img_path, files_name{i}));
        pic2 = pic2(:,:,1);
        result = dgradfunc(pic2);
    end
    result = threshold_top_percent(result, 0.001);
    file_str = fullfile(output_path, files_name{i});
    imwrite(mat2gray(result), file_str);
end

function [registered_img, translation_x, translation_y] = registration_1(I_ref, I_moving)
    % 计算刚性配准（平移 + 旋转），但去除旋转分量
    optimizer = registration.optimizer.RegularStepGradientDescent;
    optimizer.MinimumStepLength = 1e-6;
    optimizer.MaximumIterations = 30;
    metric = registration.metric.MeanSquares;

    % 估计变换矩阵
    tform = imregtform(I_moving, I_ref, 'rigid', optimizer, metric);
    
    % 去除旋转分量
    A = tform.T;
    tform.T = A;

    % 执行图像配准
    output_view = imref2d(size(I_ref));
    registered_img = imwarp(I_moving, tform, 'OutputView', output_view);

    translation_x = A(3,1);
    translation_y = A(3,2);
end

function I_thresh = threshold_top_percent(I, percent)
    I = double(I);
    sorted_vals = sort(I(:), 'descend');
    idx = round(percent * numel(sorted_vals));
    T = sorted_vals(idx);
    mask = I >= T;
    I_thresh = I .* mask;
end

function [SCR_map, TEMPORAL_map] = Morphology_analysis(filename, S_map_B, path, thr, frame_d)
% figure(); imshow(S_map_B),impixelinfo;
%%构造时序数列
    filename_org = filename;
    len = length(filename_org) - 4;
    index_filename = str2num(filename_org(1:len));
    D = zeros(256, 256, 2* thr + 1);
    for t = 1:2* thr+1
        tmp_index = t+index_filename-thr-1;
        tmp = imread(fullfile(path, string(strcat(num2str(tmp_index), '.bmp'))));
        D(:,:,t) = tmp(:,:,1);
    end
    S_map = D(:,:,thr + 1);

[row,col]=size(S_map_B);
% 
Conn_struct = bwconncomp(S_map_B); 
% 
SCR_map = zeros(row,col);
TEMPORAL_map = zeros(row,col);

%% 
for i=1:Conn_struct.NumObjects
    Conn_temp = Conn_struct.PixelIdxList(i);
    Conn_index = [Conn_temp{1,1}];          % 
    [Conn_num, no_use]= size(Conn_index);    % 
    clear Conn_position;                    % 
    
    Conn_position = zeros(Conn_num, 2);
	for j=1:Conn_num                        
        Conn_position(j,1) = mod(Conn_index(j,1),row);                      % 
        if Conn_position(j,1)==0
            Conn_position(j,1) = row;
        end
        Conn_position(j,2) = ((Conn_index(j,1)-Conn_position(j,1))/row) +1; %            
    end
    scale_factor = 5;

    win_size = scale_factor*3;             % 
    r_win = (win_size-1)/2;
    

    %% 局部对比度图生成
    for j=1:Conn_num 
        %
        img_tmp = frame_d;
        if ( (Conn_position(j,1)-r_win)>0 & (Conn_position(j,1)+r_win)<=row & (Conn_position(j,2)-r_win)>0 & (Conn_position(j,2)+r_win)<=col) ==1
            patch_SCR = img_tmp(Conn_position(j,1)-r_win:Conn_position(j,1)+r_win,Conn_position(j,2)-r_win:Conn_position(j,2)+r_win); 
            result = dgradfunc(patch_SCR);
            SCR_map(Conn_position(j,1),Conn_position(j,2)) = max(result(:));
        end
    end

    %% 时域对比度图生成
    r_win_tmporal = 2;
    
    channels = size(D, 3);
    list_max = zeros(channels, 1);
    for j = 1:Conn_num
        if ( (Conn_position(j,1)-r_win)>0 & (Conn_position(j,1)+r_win)<=row & (Conn_position(j,2)-r_win)>0 & (Conn_position(j,2)+r_win)<=col) ==1
             for t = 1:channels   
                img = D(:,:,t);
                patch_TEMPORAL = img(Conn_position(j,1)-r_win_tmporal:Conn_position(j,1)+r_win_tmporal,Conn_position(j,2)-r_win_tmporal:Conn_position(j,2)+r_win_tmporal);
                max_tep = mean(patch_TEMPORAL(:));
                list_max(t) = max_tep;
             end
             TEMPORAL_map(Conn_position(j,1),Conn_position(j,2)) = (max(list_max) - min(list_max))^2 * S_map_B(Conn_position(j,1),Conn_position(j,2));

        end    
    end

end

% SCR_map = mat2gray(SCR_map);
% TEMPORAL_map = mat2gray(TEMPORAL_map);
end

function re = dgradfunc(img)
[row, col] = size(img);
len = 3;
nr = 5;
nc = 5;
leny = len*nr;
lenx = len*nc;
op = zeros(leny, lenx, nr*nc);
for ii = 1:nr*nc
temp = zeros(len*nr, len*nc);
[r, c] = ind2sub([nr, nc], ii);
if ii ==13
sigma = 0.1772; % 控制高斯函数的标准差
gaussianFilter = fspecial('gaussian', [len, len], sigma);
temp((r-1)*len + 1:r*len, (c-1)*len + 1:c*len) = gaussianFilter;
else
temp((r-1)*len + 1:r*len, (c-1)*len+1:c*len) = 1;
end
temp = temp/sum(temp(:));
temp = temp';
op(:, :, ii) = temp;
end
gimg = zeros(row, col, nr*nc);
for ii = 1:nr*nc
gimg(:, :, ii) = imfilter(img, op(:,:, ii), 'replicate');
end
m2 = gimg;
m2(:,:, [1:6, 10:11, 13, 15:16, 20:25]) = [];
t2 = repmat(gimg(:, :, 13), 1,1,8) - m2;
t2(t2<=0) = 0;
d2 = min(t2, [], 3); %%% 内层的已经计算出来了
m1 = gimg;
m1(:,:, [7:9, 12:14, 17:19]) = [];
t1 = abs(repmat(gimg(:, :, 13), 1,1,16) - m1);
d1 = mean(t1, 3);
d1 = d1./(gimg(:, :, 13)+0.1);
re = (d1).*d2;
end

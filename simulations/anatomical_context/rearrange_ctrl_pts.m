function [ctrl_pts2,ctrl_pts] = rearrange_ctrl_pts(ctrl_pts)

ctrl_pts(:,:,1) = flipdim(ctrl_pts(:,:,1),1);
ctrl_pts(:,:,2) = flipdim(ctrl_pts(:,:,2),1);
ctrl_pts(:,:,3) = flipdim(ctrl_pts(:,:,3),1);
ctrl_pts(:,:,14) = flipdim(ctrl_pts(:,:,14),1);


ctrl_pts(1,:,1) = ctrl_pts(1,:,2);
ctrl_pts(1,:,3) = ctrl_pts(1,:,2);
ctrl_pts(1,:,4) = ctrl_pts(1,:,2);
ctrl_pts(1,:,5) = ctrl_pts(1,:,2);

ctrl_pts(1,:,6) = ctrl_pts(1,:,7);
ctrl_pts(1,:,9) = ctrl_pts(1,:,7);

ctrl_pts(1,:,8) = ctrl_pts(1,:,10);
ctrl_pts(1,:,11) = ctrl_pts(1,:,10);
ctrl_pts(1,:,12) = ctrl_pts(1,:,10);

ctrl_pts(1,:,13) = [295,390];
ctrl_pts(1,:,14) = ctrl_pts(1,:,13);




ctrl_pts(:,:,1) = ctrl_pts(:,:,2);
ctrl_pts(:,:,3) = ctrl_pts(:,:,2);
ctrl_pts(:,:,4) = ctrl_pts(:,:,5);
ctrl_pts(:,:,6) = ctrl_pts(:,:,7);
ctrl_pts(:,:,8) = ctrl_pts(:,:,10);
ctrl_pts(:,:,14) = ctrl_pts(:,:,13);

ctrl_pts(end,:,8) = ctrl_pts(end,:,9);

ctrl_pts(:,:,14) = ctrl_pts(:,:,11);
ctrl_pts(1,:,end) = ctrl_pts(1,:,13);

ctrl_pts2(:,:,1) = ctrl_pts(:,:,1); ctrl_pts2(:,:,2) = ctrl_pts(:,:,4);
ctrl_pts2(:,:,3) = ctrl_pts(:,:,6); ctrl_pts2(:,:,4) = ctrl_pts(:,:,9);
ctrl_pts2(:,:,5) = ctrl_pts(:,:,10); ctrl_pts2(:,:,6) = ctrl_pts(:,:,11);
ctrl_pts2(:,:,7) = ctrl_pts(:,:,14); ctrl_pts2(:,:,8) = ctrl_pts(:,:,13);
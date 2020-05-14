function [txx,txy,txz,tyy,tyz,tzz]=sdr2mt(strike,dip,rake,mo)
    % Function to convert strike, dip, rake and scalar moment to the 6 MT
    % components
    txx=-mo*(sind(dip)*cosd(rake)*sind(2*strike)+sind(2*dip)*sind(rake)*(sind(strike))^2);
    txy=mo*(sind(dip)*cosd(rake)*cosd(2*strike)+0.5*sind(2*dip)*sind(rake)*sind(2*strike));
    txz=-mo*(cosd(dip)*cosd(rake)*cosd(strike)+cosd(2*dip)*sind(rake)*sind(strike));
    tyy=mo*(sind(dip)*cosd(rake)*sind(2*strike)-sind(2*dip)*sind(rake)*(cosd(strike)^2));
    tyz=-mo*(cosd(dip)*cosd(rake)*sind(strike)-cosd(2*dip)*sind(rake)*cosd(strike));
    tzz=mo*(sind(2*dip)*sind(rake));
end

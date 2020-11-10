function [H] = get_report_heading(order,my_title)

import mlreportgen.dom.*;
import mlreportgen.report.*;
iptgetpref('ImshowInitialMagnification');

switch order
    
    case 1
        
        % Report Heading 1 
        H = Heading(1,my_title);
        H.FontSize="32pt";
        H.Color='#00b050';
        H.FontFamilyName='Arial';

    case 2
        
        % Report Heading 2
        H = Heading(2,my_title);
        H.FontSize="24pt";
        H.Color='#809ec2';
        H.FontFamilyName='Arial';
        
    case 3

        % Report Heading 3
        H = Heading(3,my_title);
        H.FontSize = "14pt";
        H.Color = '#404040';
        H.FontFamilyName = 'Arial';
        
    case 4

        % Report Heading 4
        H = Heading(4,my_title);
        H.FontSize = "10pt";
        H.Color = '#404040';
        H.FontFamilyName = 'Arial';
        H.Bold = false;
        
end

end
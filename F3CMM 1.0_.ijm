// Calcium calibration using FLIM data
// Using Fiji, bioformats need to be supported
//
// 2014-01-31
// Romain Laine https://github.com/Romain-Laine

Border = 0.4;               // in fraction of the width
ShiftThresh = 0.01;         // for the estimation of the IRF shift 0.2 = 20% of max
Apply_corr_factor = true;
Apply_CMF3fusion = true;

ApplyCalibration = false;

setBatchMode(true);
run("Bio-Formats Macro Extensions");

close("*FLIM*");
close("*stack - raw data*");
close("*stack - RGB*");
print("\\Clear");
run("Clear Results");

GUI_name = "F3CMM 1.0";

print(GUI_name);
print("------------------------------------------------------------------------------------------------");

getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
if (minute < 10) print("Date: "+year+"-"+month+"-"+dayOfMonth+" Time: "+hour+":0"+minute);
else print("Date: "+year+"-"+month+"-"+dayOfMonth+" Time: "+hour+":"+minute);


print("------------------------------------------------------------------------------------------------");
print("Developed by Dr. Romain Laine from the Laser Analytics group, University of Cambridge");
print("https://github.com/Romain-Laine");
print("Twitter handle: @LaineBioImaging");

Extension_type = newArray(2);
Extension_type[0] = ".sdt";
Extension_type[1] = ".tif";

Sensor_type = newArray(2);
Sensor_type[0] = "FRET sensor";
Sensor_type[1] = "Direct lifetime sensor";

Smoothing_items = newArray(5);
Smoothing_items[0] = "None";
Smoothing_items[1] = "1 (3x3 pixels)";
Smoothing_items[2] = "2 (5x5 pixels)";
Smoothing_items[3] = "3 (7x7 pixels)";
Smoothing_items[4] = "4 (9x9 pixels)";

Binning_items = newArray(6);
Binning_items[0] = "None";
Binning_items[1] = "2";
Binning_items[2] = "4";
Binning_items[3] = "8";
Binning_items[4] = "16";
Binning_items[5] = "Global binning";

// Choose the calibration method and parameters
IRF_type = newArray(3);
IRF_type[0] = "Scatterer";
IRF_type[1] = "Reference";
IRF_type[2] = "None";

// Create the dialog window
Dialog.create(GUI_name+" - Dataset");
Dialog.addMessage("----- Fluorescence lifetime dataset -----");
Dialog.addChoice("Data file extension ", Extension_type, Extension_type[0]);
Dialog.addString("Data filename contains", "Data");
Dialog.addChoice("Smoothing ", Smoothing_items, Smoothing_items[1]);
Dialog.addChoice("Binning ", Binning_items, Binning_items[0]);
Dialog.addNumber("Threshold parameter ", 1, 2, 5, "");
Dialog.addMessage("----- Instrument response function dataset -----");
Dialog.addChoice("IRF file extension ", Extension_type, Extension_type[0]);
Dialog.addString("Data filename contains", "IRF");
Dialog.addChoice("IRF type", IRF_type);
Dialog.show();

DataExtension = Dialog.getChoice();
Data_name = Dialog.getString();
Smoothing_Choice = Dialog.getChoice();
Binning_Choice = Dialog.getChoice();
Thresh_param = Dialog.getNumber();
IRFExtension = Dialog.getChoice();
IRF_name = Dialog.getString();
IRFtype_choice = Dialog.getChoice();
print("IRF type: "+IRFtype_choice);

if (Binning_Choice == "Global binning") Smoothing_Choice = Smoothing_items[0];  // in case of glabal binning, set automatically the smoothing to None


print("--------------------------------");
print("Choice of Data processing");
print("Smoothing: "+Smoothing_Choice);
print("Binning: "+Binning_Choice);
print("Threshold parameter: "+Thresh_param);
print("IRF type: "+IRFtype_choice);

// Choose the correct directory
Data_dir = getDirectory("Choose Source Directory :");

// Check that the files are all there
print("--------------------------------");
print("Source directory:");
print(Data_dir);

// Get the data file names
Stack_filename_list = Make_data_file_list(Data_dir,Data_name, DataExtension);
print("Number of files containing *"+Data_name+"*"+DataExtension+": "+lengthOf(Stack_filename_list));

// Get the information from the first FLIM dataset selected (not the IRF as it's not always used)
print("--------------------------------");
print("Getting Metadata...");
OpenFLIMdata(Stack_filename_list[0]);
rename("Dataset for Metada - stack");
getDimensions(width, height, channels, slices, frames);
print("Number of bins: "+frames);
t_window = GetTimeWindow(Stack_filename_list[0]);
print("Time window obtained from metadata: "+t_window+" ns");

selectWindow("Dataset for Metada - stack");
run("Z Project...", "start=1 stop="+frames+" projection=[Sum Slices]");
rename("Dataset for Metada - Total intensity");

CalculateMask("Dataset for Metada - Total intensity", Thresh_param);
imageCalculator("Multiply create 32-bit stack", "Dataset for Metada - stack","Dataset for Metada - Total intensity mask");
rename("Dataset for Metada - stack merged");

Metadata_decay = GetDecayData("Dataset for Metada - stack merged", false, 0, 0);
//Array.show(Metadata_decay);
BGbin_max = GetBackgroundPosition(Metadata_decay);

close("Dataset for Metada - Total intensity");
close("Dataset for Metada - stack");
close("Dataset for Metada - stack merged");
close("Dataset for Metada - Total intensity mask");

t = newArray(frames);
for (i=0; i<frames; i++) t[i] = i*(t_window/frames);

// Get the IRF file name and its metadata
if (IRFtype_choice != "None"){
	IRF_filename_list = Make_data_file_list(Data_dir,IRF_name, IRFExtension);
	if (lengthOf(IRF_filename_list) > 1) exit("More than one file containing *"+IRF_name+"*"+IRFExtension+" were found.");
	IRF_filename = IRF_filename_list[0];
	print("----------------------------------------------");
	print("Opening IRF file...");
	OpenFLIMdata(IRF_filename);
	rename("FLIM IRF stack");
	getDimensions(width, height, channels, slices, frames_IRF);
	t_window_IRF = GetTimeWindow(IRF_filename);
	print("Time window obtained from IRF metadata: "+t_window_IRF+" ns");
	IRF_decay = GetDecayData("FLIM IRF stack", false, 0, 0);
	BGbin_max_IRF = GetBackgroundPosition(IRF_decay);
	t_IRF = newArray(frames_IRF);
	for (i=0; i<frames_IRF; i++) t_IRF[i] = i*(t_window_IRF/frames_IRF);}
else {
	t_max = FindPositionOfMaximum(t, Metadata_decay);
	print("Position of maximum in decay: "+t_max+" ns");}

// --------------------------------------------------------------------------------------------------------------------------------------
// Create the dialog window
Dialog.create(GUI_name+" - Background correction");
Dialog.addMessage("----- Fluorescence lifetime dataset -----");
Dialog.addCheckbox("Remove decay background",true);
Dialog.addNumber("From bin #", 1, 0, 2, "");
Dialog.addNumber("to bin #", BGbin_max, 0, 2, "");

if (IRFtype_choice != "None"){
	Dialog.addMessage("----- Instrument response function dataset -----");
	Dialog.addCheckbox("Remove decay background",true);
	Dialog.addNumber("From bin #", 1, 0, 2, "");
	Dialog.addNumber("to bin #", BGbin_max_IRF, 0, 2, "");}
	
Dialog.show();

Remove_bg = Dialog.getCheckbox();
BGbin_min = Dialog.getNumber();
BGbin_max = Dialog.getNumber();

if (IRFtype_choice != "None"){
	Remove_bg_IRF = Dialog.getCheckbox();
	BGbin_min_IRF = Dialog.getNumber();
	BGbin_max_IRF = Dialog.getNumber();}

		
print("----------------------------------------------");
print("- Dataset background:");
if (Remove_bg == false) print("Background removal: None");
if (Remove_bg == true)  print("Background removal: Yes (calculated from the average of bin #"+BGbin_min+" to bin #"+BGbin_max+")");

if (IRFtype_choice != "None"){
	print("- IRF background:");
	if (Remove_bg_IRF == false) print("Background removal: None");
	if (Remove_bg_IRF == true)  {
		print("Background removal: Yes (calculated from the average of bin #"+BGbin_min_IRF+" to bin #"+BGbin_max_IRF+")");
		IRF_decay = GetDecayData("FLIM IRF stack", Remove_bg_IRF, BGbin_min_IRF, BGbin_max_IRF);}
}


// --------------------------------------------------------------------------------------------------------------------------------------
IRF_shift_Methods = newArray(2);
IRF_shift_Methods[0] = "From dataset";
IRF_shift_Methods[1] = "User set";

DefaultName = File.getName(Stack_filename_list[0]);
DefaultName = substring(DefaultName, 0, lengthOf(DefaultName)-lengthOf(DataExtension));

if (IRFtype_choice == "Scatterer"){
	Dialog.create(GUI_name+" - IRF shift correction");
	Dialog.addMessage("----- IRF (Scatterer) settings -----");
	Dialog.addCheckbox("Estimate IRF shift",true);
	Dialog.addChoice("IRF shift estimation method", IRF_shift_Methods);
	Dialog.addMessage("- For calibration using a dataset:");
	Dialog.addString("Data contains", DefaultName);
	Dialog.addMessage("- For calibration using a fixed shift:");
	Dialog.addNumber("IRF shift", 0.2, 2, 5, "ns");
	Dialog.show();

	Estimate_IRF_shift = Dialog.getCheckbox();
	IRF_shift_Method_choice = Dialog.getChoice();
	IRF_shiftDataset_name = Dialog.getString();
	Fixed_IRF_shift = Dialog.getNumber();
	print("----------------------------------------------");
	if (Estimate_IRF_shift == true) print("IRF shift correction: Yes");
	if (Estimate_IRF_shift == false) print("IRF shift correction: No");}
	
if (IRFtype_choice == "Reference"){
	Dialog.create(GUI_name+" - IRF settings");
	Dialog.addMessage("----- IRF (Reference) settings -----");
	Dialog.addNumber("Reference lifetime", 0, 3, 5, "ns");
	Dialog.addCheckbox("Estimate IRF shift",true);
	Dialog.addChoice("IRF shift estimation method", IRF_shift_Methods);
	Dialog.addMessage("- For calibration using a dataset:");
	Dialog.addString("Data contains", DefaultName);
	Dialog.addMessage("- For calibration using a fixed shift:");
	Dialog.addNumber("IRF shift", 0.2, 2, 5, "ns");
	Dialog.show();

	Ref_lifetime = Dialog.getNumber();
	Estimate_IRF_shift = Dialog.getCheckbox();
	IRF_shift_Method_choice = Dialog.getChoice();
	IRF_shiftDataset_name = Dialog.getString();
	Fixed_IRF_shift = Dialog.getNumber();}


if (IRFtype_choice == "None") Estimate_IRF_shift = false;
else {
	// Check for the dataset used to estimate the IRF shift
	if ((Estimate_IRF_shift == true) && (IRF_shift_Method_choice == "From dataset")){
		IRFshift_filename_list = Make_data_file_list(Data_dir,IRF_shiftDataset_name, DataExtension);
		if (lengthOf(IRFshift_filename_list) > 1) exit("More than one file containing *"+IRF_shiftDataset_name+"*"+DataExtension+" were found.");
		IRFshift_filename = IRFshift_filename_list[0];}
}

// --------------------------------------------------------------------------------------------------------------------------------------

// Choose the calibration method and parameters
Methods = newArray(2);
Methods[0] = "Average lifetime";
Methods[1] = "Fractions";

Calib_method = newArray(2);
Calib_method[0] = "Tau_min/Tau_max";
Calib_method[1] = "high/low [Ca2+] dataset";

if (ApplyCalibration == true){
// Create the dialog window
Dialog.create(GUI_name+" - Calibration");
Dialog.addMessage("----- Calibration parameters -----");
Dialog.addChoice("Sensor type ", Sensor_type, Sensor_type[0]);
if (IRFtype_choice != "None") Dialog.addChoice("Analysis method", Methods);
Dialog.addNumber("Dissociation constant", 170.0, 2, 5, "uM");
Dialog.addNumber("Hill's coefficient", 1, 2, 5,"");

Dialog.addMessage("----- Zero [Ca2+] and saturating [Ca2+] datapoints -----");
Dialog.addChoice("Zero/saturating [Ca2+] method", Calib_method);
Dialog.addMessage("- For calibration using Tau_min and Tau_max:");
Dialog.addNumber("Tau_min", 1.012, 3, 5, "ns");
Dialog.addNumber("Tau_max", 2.804, 3, 5, "ns");
Dialog.addMessage("- For calibration using high and low [Ca2+] dataset:");
Dialog.addString("Saturating [Ca2+] data image contains", "highCa");
Dialog.addString("No [Ca2+] data image contains", "lowCa");
Dialog.show();

Sensor_type_choice = Dialog.getChoice();
if (IRFtype_choice != "None") AnalysisMethod_choice = Dialog.getChoice();
kD = Dialog.getNumber();
h = Dialog.getNumber();
CalibrationMethod_choice = Dialog.getChoice();
Tau_min = Dialog.getNumber();
Tau_max = Dialog.getNumber();
HighCa_name = Dialog.getString();
LowCa_name = Dialog.getString();
AnalysisMethod_choice = "Average lifetime"; //   ----------------- WARNING! only uses the average lifetime currently

if (CalibrationMethod_choice == "high/low [Ca2+] dataset"){
	HighCa_filename_list = Make_data_file_list(Data_dir, HighCa_name, DataExtension);
	if (lengthOf(HighCa_filename_list) > 1) exit("More than one file containing *"+HighCa_name+"*"+DataExtension+" were found.");
	HighCa_filename = HighCa_filename_list[0];

	LowCa_filename_list = Make_data_file_list(Data_dir, LowCa_name, DataExtension);
	if (lengthOf(LowCa_filename_list) > 1) exit("More than one file containing *"+LowCa_name+"*"+DataExtension+" were found.");
	LowCa_filename = LowCa_filename_list[0];

	if (IRFtype_choice != "None"){
		IRF_Calib_filename_list = Make_data_file_list(Data_dir, IRF_name_min_max, IRFExtension);
		if (lengthOf(IRF_Calib_filename_list) > 1) exit("More than one file containing *"+IRF_name_min_max+"*"+IRFExtension+" were found.");
		IRF_Calib_filename = IRF_Calib_filename_list[0];}
}

} // if (ApplyCalibration == true)
else {
	CalibrationMethod_choice = "None";}

// Create the dialog window for output
Dialog.create(GUI_name+" - Output");
Dialog.addMessage("----- Display output -----");
Dialog.addCheckbox("Display decays and masks",false);
Dialog.addCheckbox("Display lifetime image",true);
Dialog.addCheckbox("Display intensity-merged lifetime image",true);

if (ApplyCalibration == true){
	Dialog.addCheckbox("Display [Ca2+] image",false);
	Dialog.addCheckbox("Display intensity-merged [Ca2+] image",false);}

Dialog.addCheckbox("Display total intensity image",true);

Dialog.addMessage("----- Lifetime image scale -----");
Dialog.addCheckbox("Auto",false);
Dialog.addNumber("Minimum", 0, 1, 2, " ns");
Dialog.addNumber("Maximum", 5, 1, 2, " ns");

if (ApplyCalibration == true){
	Dialog.addMessage("----- [Ca2+] image scale bar -----");
	Dialog.addCheckbox("Auto",false);
	Dialog.addNumber("Minimum", 50, 1, 7, " uM");
	Dialog.addNumber("Maximum", 5000, 1, 7, " uM");
	Dialog.addCheckbox("Display [Ca2+] in log scale",false);
	Dialog.addCheckbox("Remove outliers",true);}

Dialog.addMessage("----- Saving output -----");
Dialog.addCheckbox("Save single images",true);
Dialog.addCheckbox("Save stacks",true);
Dialog.addString("Output folder name", "Results");
Dialog.addCheckbox("Large batch mode",true);
Dialog.show();

Keep_nonResults = Dialog.getCheckbox();
Disp_TauImage = Dialog.getCheckbox();
Disp_TauImage_merged = Dialog.getCheckbox();

if (ApplyCalibration == true){
	Disp_CaImage = Dialog.getCheckbox();
	Disp_CaImage_merged = Dialog.getCheckbox();}
else {
	Disp_CaImage = false;
	Disp_CaImage_merged = false;}

Disp_IntImage = Dialog.getCheckbox();

Scale_Tau_Auto = Dialog.getCheckbox();
Scale_Tau_min = Dialog.getNumber();
Scale_Tau_max = Dialog.getNumber();

if (ApplyCalibration == true){
	Scale_Ca_Auto = Dialog.getCheckbox();
	Scale_Ca_min = Dialog.getNumber();
	Scale_Ca_max = Dialog.getNumber();
	Use_log_scale = Dialog.getCheckbox();
	Remove_outliers = Dialog.getCheckbox();}


Save_SingleImages_on = Dialog.getCheckbox();
Save_Stacks_on = Dialog.getCheckbox();

Output_folder_name = Dialog.getString();
LargeBatchMode = Dialog.getCheckbox();

if (LargeBatchMode == 1) Keep_nonResults = 0;

print("--------------------------------");
if (Save_SingleImages_on == 0) Save_yn = "No";
if (Save_SingleImages_on == 1) Save_yn = "Yes";
print("Save single images selected? "+Save_yn);

if (Save_Stacks_on == 0) Save_yn = "No";
if (Save_Stacks_on == 1) Save_yn = "Yes";
print("Save stacks selected? "+Save_yn);

// --------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------

if (IRFtype_choice == "Scatterer"){
	// Crop the IRF stack to the 1 pixel in the corner that contains the information
	selectWindow("FLIM IRF stack");
	makeRectangle(0, 0, 1, 1);
	run("Crop");
	Ref_lifetime = 0; // also set the reference lifetime to zeros
}
	
if (IRFtype_choice == "Reference"){
	// Apply global binning to the image if it's a reference image
	selectWindow("FLIM IRF stack");
	run("Bin...", "x=&width y=&height z=1 bin=Sum");}

// Calculate the IRF shift
if (Estimate_IRF_shift == true){
	if (IRF_shift_Method_choice == "From dataset"){
		print("--------------------------------");
		print("Estimating the IRF shift...");
			
		// Estimate the T_thresh for the IRF shift estimation 
		T_thresh_IRF = GetTimeThresh(t_IRF, IRF_decay, ShiftThresh);
		T_thresh_bin_IRF = TimeThresh2Bin_position(T_thresh_IRF, t_window_IRF, frames_IRF);

		if (IRFtype_choice == "Scatterer")	Shift_correction = Shift_correction_factor(t_IRF, IRF_decay);
		else Shift_correction = 0;

		// Load IRF shift estimation dataset 
		print("Opening dataset for IRF shift estimation..."); 
		OpenFLIMdata(IRFshift_filename);
		rename("IRF shift stack");

		getDimensions(width, height, channels, slices, frames);
		run("Z Project...", "start=1 stop="+frames+" projection=[Sum Slices]");
		rename("FLIM IRF shift dataset - Total intensity");
	
   		CalculateMask("FLIM IRF shift dataset - Total intensity",Thresh_param);
   		//CalculateMask2("FLIM High Ca stack");
		imageCalculator("Multiply create 32-bit stack", "IRF shift stack","FLIM IRF shift dataset - Total intensity mask");
		rename("IRF shift stack merged");
		IRFshift_decay = GetDecayData("IRF shift stack merged", Remove_bg, BGbin_min, BGbin_max);
		T_thresh_IRFshift = GetTimeThresh(t, IRFshift_decay, ShiftThresh);
		T_thresh_bin = TimeThresh2Bin_position(T_thresh_IRFshift, t_window, frames);
		Shift = T_thresh_IRFshift - T_thresh_IRF - Shift_correction;
		
		if (isNaN(Shift)) Shift = 0;
		print("IRF shift estimate: "+d2s(Shift,3)+" ns");

		close("FLIM IRF shift dataset - Total intensity");
		close("IRF shift stack");
		close("IRF shift stack merged");
		close("IRF shift stack mask");}

	else if (IRF_shift_Method_choice == "3User set") {
		Shift = Fixed_IRF_shift;
		print("IRF shift was set by user to "+Fixed_IRF_shift+" ns.");}
} // end of if (Estimate_IRF_shift == true) 

// Alternaively shift should be set to zero if it was set not to be estimated 
if (Estimate_IRF_shift == false) Shift = 0;

if (IRFtype_choice != "None"){
	// Calculate the lifetime image of the one pixel image
	// Setting the lifetime of the Tau_irf in the following function is a problem for reference lifetime used as IRF with long lifetime. 
	CalculateAvLifetime("FLIM IRF stack",t_window_IRF, Remove_bg_IRF, BGbin_min_IRF, BGbin_max_IRF, 0, false, false);
	rename("IRF lifetime image");

	Tau_IRF = getPixel(0,0) - Ref_lifetime + Shift;
	print("IRF position: "+d2s(Tau_IRF,3)+" ns");

	// Adaptive window ___BOTCH BOTCH BOTCH !!!!
	//Tau_IRF = 3.353;

	close("IRF lifetime image");
	close("FLIM IRF stack");
	close("FLIM IRF stack_original");
} // end of IRF choice != "None"
else Tau_IRF = t_max;  // set Tau_IRF to zero if there is no IRF ("None")


// Open the Tau_min and Tau_max datasets
if (CalibrationMethod_choice == "high/low [Ca2+] dataset"){

	// Is the IRF the same for the High / low calcium and for the datasets?
	if (IRFtype_choice != "None"){
		if (IRF_name_min_max == IRF_name){
			Tau_IRF_Cal = Tau_IRF;
			t_window_IRF_Calib = t_window_IRF;
			t_IRF_Calib = t_IRF;
			if (Keep_nonResults == 1) IRF_Cal_decay = IRF_decay;}
		
		else{
			print("--------------------------------");
			print("Opening IRF Calibration dataset...");
			OpenFLIMdata(IRF_Calib_filename);
			rename("FLIM Cal IRF stack");
			getDimensions(width, height, channels, slices, frames);
			t_window_IRF_Calib = GetTimeWindow(IRF_Calib_filename);

			t_IRF_Calib = newArray(frames);
			for (i=0;i<frames;i++) t_IRF_Calib[i] = i*(t_window_IRF_Calib/frames);

			selectWindow("FLIM Cal IRF stack");
			makeRectangle(0, 0, 1, 1);
			run("Crop");

			if (Keep_nonResults == 1) IRF_Cal_decay = GetDecayData("FLIM Cal IRF stack", Remove_bg_IRF, BGbin_min_IRF, BGbin_max_IRF);
			
			// Setting the lifetime of the Tau_irf in the following function is a problem for reference lifetime used as IRF with long lifetime. 
			CalculateAvLifetime("FLIM Cal IRF stack",t_window_IRF_Calib, Remove_bg_IRF, BGbin_min_IRF, BGbin_max_IRF, 0, false, false);
			rename("IRF lifetime image");
			Tau_IRF_Cal = getPixel(0,0);
			
			// Here the background and shift correction is missing !!!
			//
			//
			//			

			
			print("IRF position for calibration: "+d2s(Tau_IRF_Cal,3)+" ns");
			close("IRF lifetime image");}
	} // end of if (IRFtype_choice != "None")
	else Tau_IRF_Cal = 0;
	
	// Load High Calcium dataset and calculate Tau_min 
	print("--------------------------------");
	print("Opening High Calcium dataset...");
	OpenFLIMdata(HighCa_filename);
	rename("FLIM High Ca stack");
	getDimensions(width, height, channels, slices, frames);
	run("Z Project...", "start=1 stop="+frames+" projection=[Sum Slices]");
	rename("FLIM High Ca stack - Total intensity");

    CalculateMask("FLIM High Ca stack - Total intensity",Thresh_param);
    //CalculateMask2("FLIM High Ca stack - Total intensity");
    imageCalculator("Multiply create 32-bit stack", "FLIM High Ca stack","FLIM High Ca stack - Total intensity mask");
    rename("FLIM High Ca stack merged");

    if (Keep_nonResults == 1) HighCa_decay = GetDecayData("FLIM High Ca stack merged", Remove_bg, BGbin_min, BGbin_max);

    selectWindow("FLIM High Ca stack merged");
    getDimensions(width, height, channels, slices, frames);
    run("Bin...", "x=&width y=&height z=1 bin=Sum");
    
	CalculateAvLifetime("FLIM High Ca stack merged",t_window, Remove_bg, BGbin_min, BGbin_max, Tau_IRF_Cal, Apply_corr_factor, Apply_CMF3fusion);
	rename("High Ca lifetime image");

	if (Sensor_type_choice == "FRET sensor") Tau_min = getPixel(0,0);
	if (Sensor_type_choice == "Direct lifetime sensor") Tau_max = getPixel(0,0);

	close("High Ca lifetime image");
	close("FLIM High Ca stack");
	close("FLIM High Ca stack merged");
	close("FLIM High Ca stack merged_original");

	// Load Low Calcium dataset and calculate Tau_max
	print("--------------------------------");
	print("Opening Low Calcium dataset...");
	OpenFLIMdata(LowCa_filename);
	rename("FLIM Low Ca stack");
	getDimensions(width, height, channels, slices, frames);
	run("Z Project...", "start=1 stop="+frames+" projection=[Sum Slices]");
	rename("FLIM Low Ca stack - Total intensity");

    CalculateMask("FLIM Low Ca stack - Total intensity",Thresh_param);
    //CalculateMask2("FLIM Low Ca stack - Total intensity");
    imageCalculator("Multiply create 32-bit stack", "FLIM Low Ca stack","FLIM Low Ca stack - Total intensity mask");
    rename("FLIM Low Ca stack merged");

    if (Keep_nonResults == 1) LowCa_decay = GetDecayData("FLIM Low Ca stack merged", Remove_bg, BGbin_min, BGbin_max);

    selectWindow("FLIM Low Ca stack merged");
    getDimensions(width, height, channels, slices, frames);
    run("Bin...", "x=&width y=&height z=1 bin=Sum");
	CalculateAvLifetime("FLIM Low Ca stack merged", t_window, Remove_bg, BGbin_min, BGbin_max, Tau_IRF_Cal, Apply_corr_factor, Apply_CMF3fusion);
	rename("Low Ca lifetime image");
	
	if (Sensor_type_choice == "FRET sensor") Tau_max = getPixel(0,0);
	if (Sensor_type_choice == "Direct lifetime sensor") Tau_min = getPixel(0,0);
	
	close("Low Ca lifetime image");
	close("FLIM Low Ca stack");
	close("FLIM Low Ca stack merged");
	close("FLIM Low Ca stack merged_original");

	// This is now performed in the CalculateAvLifetime function
	// Tau_min = Tau_min - Tau_IRF_Cal;
	// Tau_max = Tau_max - Tau_IRF_Cal;
	
	} // end of 
if (CalibrationMethod_choice == Calib_method[1])

if(Disp_CaImage == 1 || Disp_CaImage_merged == 1){
	print("--------------------------------");
	print("Calibration parameters:");
	print("Analysis method: "+AnalysisMethod_choice);
	print("Calibration method: "+CalibrationMethod_choice);
	print("kD = "+kD+" uM");
	print("h = "+h);
	print("Tau_min = "+Tau_min+" ns");
	print("Tau_max = "+Tau_max+" ns");}

if (Keep_nonResults == 0){
	close("*stack*");
	close("*mask*");
}

// Exit batch mode
 
// This allows display of the images obtained so far
if (LargeBatchMode == false) setBatchMode("exit and display");

if ((Keep_nonResults == 1) && (CalibrationMethod_choice == "high/low [Ca2+] dataset")){
	Plot.create("FLIM Calibration decays", "Time (ns)", "Intensity (AU)");
	Plot.setLimits(0, t_window, -0.1, 1.1);
	
	if (IRFtype_choice != "None"){
		Plot.setColor("black");
		Plot.add("line", t_IRF_Calib, IRF_Cal_decay);
}
		
	Plot.setColor("red");
	Plot.add("line", t,LowCa_decay);
	Plot.setColor("blue");
	Plot.add("line", t,HighCa_decay);
	Plot.show();}


print("--------------------------------");
if (Scale_Tau_Auto == 0) print("Limit of lifetime scale bar (in ns): "+d2s(Scale_Tau_min,2)+" - "+d2s(Scale_Tau_max,2));

if (ApplyCalibration == true){
if (Scale_Ca_Auto == 0){
	if (Use_log_scale == 0) print("Limit of [Ca2+] scale bar (in uM): "+d2s(Scale_Ca_min,1)+" - "+d2s(Scale_Ca_max,1));
	if (Use_log_scale == 1){
		if (Scale_Ca_min < 1) Scale_Ca_min = 1;
			LogScale_Ca_min = log(Scale_Ca_min) / log(10);
			LogScale_Ca_max = log(Scale_Ca_max) / log(10);
			print("Limit of [Ca2+] scale bar (in log10 of uM): "+LogScale_Ca_min+" - "+LogScale_Ca_max);
}
}
}


if (Save_SingleImages_on == 1 || Save_Stacks_on == 1){
// Creating the variable and directory for the output directory
	Output_dir = Data_dir + Output_folder_name + File.separator;
	File.makeDirectory(Output_dir);}
if (Save_SingleImages_on == 1) {
	Output_dir_SingleImages = Output_dir + "Single images" + File.separator;
	File.makeDirectory(Output_dir_SingleImages);}

// Bin the dataset
if (Binning_Choice == Binning_items[0]) Data_binning = 1;
if (Binning_Choice == Binning_items[1]) Data_binning = 2;
if (Binning_Choice == Binning_items[2]) Data_binning = 4;
if (Binning_Choice == Binning_items[3]) Data_binning = 8;
if (Binning_Choice == Binning_items[4]) Data_binning = 16;
// Do not apply the global binning at that stage
if (Binning_Choice == Binning_items[5]) Data_binning = 1;  
// in case of glabal binning, the smoothing was previously automatically set to None

// Select the correct kernel
if (Smoothing_Choice == Smoothing_items[1]) Kernel = "[1 1 1\n1 1 1\n1 1 1\n]";
if (Smoothing_Choice == Smoothing_items[2]) Kernel = "[1 1 1 1 1\n1 1 1 1 1\n1 1 1 1 1\n1 1 1 1 1 \n1 1 1 1 1 \n]";
if (Smoothing_Choice == Smoothing_items[3]) Kernel = "[1 1 1 1 1 1 1\n1 1 1 1 1 1 1\n1 1 1 1 1 1 1\n1 1 1 1 1 1 1\n1 1 1 1 1 1 1\n1 1 1 1 1 1 1\n1 1 1 1 1 1 1\n]";
if (Smoothing_Choice == Smoothing_items[4]) Kernel = "[1 1 1 1 1 1 1 1 1\n1 1 1 1 1 1 1 1 1\n1 1 1 1 1 1 1 1 1\n1 1 1 1 1 1 1 1 1\n1 1 1 1 1 1 1 1 1\n1 1 1 1 1 1 1 1 1\n1 1 1 1 1 1 1 1 1\n1 1 1 1 1 1 1 1 1\n1 1 1 1 1 1 1 1 1\n]";


// ------------------------------------------------------------------------------------------------------------------------------
run("Set Measurements...", "area mean standard modal min median skewness redirect=None decimal=3");

// Here we go !
print("--------------------------------");
t0 = getTime();

// Starts batch loop
// --------------------------------------------------------------------------------------------------------------------------------------
for (j=0;j<Stack_filename_list.length;j++){

t_start = getTime();

if (LargeBatchMode == false) setBatchMode(true);
  // set it back to true if that was set to "exit and display" cf line 460
print("------------------------------------------------------------------------------------------------");
File_number = j+1;

// Open dataset file 1 by 1
//t1 = getTime();
print("Opening file number: "+File_number+"/"+Stack_filename_list.length);
OpenFLIMdata(Stack_filename_list[j]);
rename("FLIM stack");
getDimensions(width, height, channels, slices, frames);
//t2 = getTime();
//dt = (t2-t1)/1000;
//print("Time elapsed openSDTimage: "+dt+" s");

// Obtain the filename without the directory appended 
StartIndex = lastIndexOf(Stack_filename_list[j], File.separator);
EndIndex = indexOf(Stack_filename_list[j], DataExtension);
Short_filename = substring(Stack_filename_list[j], StartIndex+1, EndIndex);

//t1 = getTime();
if ((Binning_Choice != Binning_items[0]) && (Binning_Choice != Binning_items[5])){
	// Apply the appropriate binning
	selectWindow("FLIM stack");
	run("Bin...", "x=&Data_binning y=&Data_binning z=1 bin=Sum");}

selectWindow("FLIM stack");
getDimensions(width, height, channels, slices, frames);

// Calculate the total intensity image before the smoothing (arguable but looks better)
selectWindow("FLIM stack");
run("Z Project...", "start=1 stop="+frames+" projection=[Sum Slices]");
rename("FLIM total intensity");

//t2 = getTime();
//dt = (t2-t1)/1000;
//print("Time elapsed binning and summing: "+dt+" s");

//t1 = getTime();
// Calculate the mask
CalculateMask("FLIM total intensity",Thresh_param);
selectWindow("FLIM total intensity mask");
rename("FLIM mask");

//t2 = getTime();
//dt = (t2-t1)/1000;
//print("Time elapsed calculate mask: "+dt+" s");

//t1 = getTime();
// Smooth the dataset using square kernel
if (Smoothing_Choice != Smoothing_items[0]){
	// Apply the kernel by convolution
	selectWindow("FLIM stack");
	run("Convolve...", "text1="+Kernel+" stack");}
//t2 = getTime();
//dt = (t2-t1)/1000;
//print("Time elapsed smooth: "+dt+" s");	

//t1 = getTime();
// Apply the mask on the stack to get rid of the dark pixels in the decay calculation
imageCalculator("Multiply create 32-bit stack", "FLIM stack","FLIM mask");
rename("FLIM stack masked");
//t2 = getTime();
//dt = (t2-t1)/1000;
//print("Time elapsed apply mask: "+dt+" s");

//t1 = getTime();
// Calculate the decay
if (Keep_nonResults == 1) Data_decay = GetDecayData("FLIM stack masked", Remove_bg, BGbin_min, BGbin_max);

// If it's not global binning
if (Binning_Choice != Binning_items[5]){
	// Calculate the lifetime image
	CalculateAvLifetime("FLIM stack masked",t_window, Remove_bg, BGbin_min, BGbin_max, Tau_IRF, Apply_corr_factor, Apply_CMF3fusion);
	rename("Average lifetime image");
}

// In the case of global binning
if (Binning_Choice == Binning_items[5]){
	selectWindow("FLIM stack masked");
    getDimensions(width, height, channels, slices, frames);
    run("Bin...", "x=&width y=&height z=1 bin=Sum");
	CalculateAvLifetime("FLIM stack masked",t_window, Remove_bg, BGbin_min, BGbin_max, Tau_IRF, Apply_corr_factor, Apply_CMF3fusion);
	rename("Single pixel lifetime image");
	Tau_single_pixel = getPixel(0,0);
	newImage("Average lifetime image", "32-bit black", width, height, 1);
	
	if (isNaN(Tau_single_pixel) == true) run("Add...", "value=NaN"); // in the case where there are no pixels in the mask
	else{
	run("Add...", "value=&Tau_single_pixel");}
	
	close("Single pixel lifetime image");
	}

//t2 = getTime();
//dt = (t2-t1)/1000;
//print("Time elapsed calculate lifetime: "+dt+" s");


//t1 = getTime();
// Now that the lifetime image has been generated, carry on
//selectWindow("Data lifetime image");
// run("Subtract...", "value="+Tau_IRF);
//rename("Average lifetime image");
//GetRidNegative("Average lifetime image");
close("FLIM stack masked");
//t2 = getTime();
//dt = (t2-t1)/1000;
//print("Time elapsed get rid of negative: "+dt+" s");

// Apply the mask

//t1 = getTime();
imageCalculator("Multiply 32-bit", "Average lifetime image","FLIM mask");
rename("Average lifetime image - Masked");
// Multiply and divide by 0 to create NaN values, this way the values outside of the mask are excluded from the histogram and statistics measurements 
imageCalculator("Divide 32-bit", "Average lifetime image - Masked","FLIM mask");
//t2 = getTime();
//dt = (t2-t1)/1000;
//print("Time elapsed apply mask to lifetime image: "+dt+" s");

rename("FLIM lifetime image - raw data (in ns)");
run("Rainbow RGB");
//run("16_colors");
 // in case we want to use a different LUT 

//t1 = getTime();
if (Scale_Tau_Auto == 1){
	run("Enhance Contrast", "saturated=0.35");
	getMinAndMax(Scale_Tau_min, Scale_Tau_max);
	print("Limit of lifetime scale bar (in ns): "+d2s(Scale_Tau_min,2)+" - "+d2s(Scale_Tau_max,2));}


// Set the scale as defined by the user
if (Scale_Tau_Auto == 0){
	selectWindow("FLIM lifetime image - raw data (in ns)");
	setMinAndMax(Scale_Tau_min, Scale_Tau_max);
}


//t2 = getTime();
//dt = (t2-t1)/1000;
//print("Time elapsed scale image: "+dt+" s");

close("Data lifetime image");
close("Average lifetime image");
close("Average lifetime image - Masked");

//t1 = getTime();
if(Disp_TauImage == 1 || Disp_TauImage_merged == 1){
	// Measure the info of the image
	selectWindow("FLIM lifetime image - raw data (in ns)");
	run("Measure");
	setResult("Label", nResults -1 , "Lifetime (ns) - "+Short_filename);
	updateResults();

	// Create and append the lifetime stack
	if (j == 0) {
		selectWindow("FLIM lifetime image - raw data (in ns)");
		run("Duplicate...", "title=[Lifetime stack - raw data (ns)]");}
	else AppendToStack("Lifetime stack - raw data (ns)", "FLIM lifetime image - raw data (in ns)");

	// Display the FLIM image with a scale bar and merged to the intensity image
	selectWindow("FLIM lifetime image - raw data (in ns)");
	Border_size = floor(Border*height);
	//print("Border size: "+Border_size+" pixels");
	AddScaleBar("FLIM lifetime image - raw data (in ns)", Border_size);
	rename("FLIM lifetime image");
} // end of if(Disp_TauImage == 1 || Disp_TauImage_merged == 1)


if(Disp_TauImage == 1){
	// Create and append the stack 
	if (j == 0) {
		selectWindow("FLIM lifetime image");
		run("Duplicate...", "title=[Lifetime stack - RGB]");}
	else AppendToStack("Lifetime stack - RGB", "FLIM lifetime image");
}
	

if(Disp_TauImage_merged == 1){
	selectWindow("FLIM total intensity");
	run("Duplicate...", "title=[FLIM total intensity - mask]");
	getMinAndMax(min, max);
	run("Divide...", "value=max");

	new_width = width + Border_size;
	run("Canvas Size...", "width=&new_width height=&height position=Center-Left zero");
	setColor(1);
	fillRect(width+1, 0, Border_size, height);
	MergeIntensityImage("FLIM lifetime image", "FLIM total intensity - mask");
	close("FLIM total intensity - mask");
		
	// Create and append the stack 
	if (j == 0) {
		selectWindow("FLIM lifetime image merged");
		run("Duplicate...", "title=[Lifetime merged stack - RGB]");}
	else AppendToStack("Lifetime merged stack - RGB", "FLIM lifetime image merged");
}


//t2 = getTime();
//dt = (t2-t1)/1000;
//print("Time elapsed lifetime images with scale bar: "+dt+" s");


if(Disp_CaImage == 1 || Disp_CaImage_merged == 1){
	// Calculation for the calibration
	// ----------------------- Pixelwise calculation -------------------------------
	selectWindow("FLIM lifetime image - raw data (in ns)");
	run("Duplicate...", "title=[FLIM Calcium concentration - raw data (in uM)]");

	for(x=0;x<width;x++){
		for(y=0;y<height;y++){
			Tau = getPixel(x,y);
			if (Tau == 0) Ca = NaN;
			else {
				if(Tau < Tau_min) Tau = Tau_min;
				if(Tau > Tau_max) Tau = Tau_max;

				if (Sensor_type_choice == "FRET sensor") 
					{T1 = Tau_max; 
					T2 = Tau_min;}
				if (Sensor_type_choice == "Direct lifetime sensor")
					{T1 = Tau_min;
					T2 = Tau_max;}
			
				// The two methods use different formulae
				if (AnalysisMethod_choice == "Average lifetime") base = kD*(T1-Tau)/(Tau-T2);
				if (AnalysisMethod_choice == "Fractions") base = kD*(T1/T2)*(T1-Tau)/(Tau-T2);
				Ca = pow(base, 1/h);}
			setPixel(x,y,Ca);
}}

	// Post-processing
	if (Remove_outliers == 1) {
	run("Remove Outliers...", "radius=2 threshold=100 which=Dark");
	run("Remove Outliers...", "radius=2 threshold=100 which=Bright");
}

	if (Use_log_scale == 0 && Scale_Ca_Auto == 1){
	selectWindow("FLIM Calcium concentration - raw data (in uM)");
	run("Enhance Contrast", "saturated=0.35");
	getMinAndMax(Scale_Ca_min, Scale_Ca_max);
	print("Limit of [Ca2+] scale bar (in uM): "+d2s(Scale_Ca_min,1)+" - "+d2s(Scale_Ca_max,1));}

	if (Use_log_scale == 0 && Scale_Ca_Auto == 0){
	selectWindow("FLIM Calcium concentration - raw data (in uM)");
	run("Enhance Contrast", "saturated=0.35");
	setMinAndMax(Scale_Ca_min, Scale_Ca_max);}

	if (Use_log_scale == 1){
	selectWindow("FLIM Calcium concentration - raw data (in uM)");
	run("Log");
	ln10 = log(10);
	run("Divide...", "value=&ln10");
	run("Enhance Contrast", "saturated=0.35");
	
	if (Scale_Ca_Auto == 0){	
		selectWindow("FLIM Calcium concentration - raw data (in uM)");
		setMinAndMax(LogScale_Ca_min, LogScale_Ca_max);}
}

	// Measure the statistics
	selectWindow("FLIM Calcium concentration - raw data (in uM)");
	run("Measure");
	setResult("Label", nResults -1 , "[Ca2+] (uM)      - "+Short_filename);
	updateResults();

		// Create and append the concentration stack 
	if (j == 0) {
		selectWindow("FLIM Calcium concentration - raw data (in uM)");
		run("Duplicate...", "title=[Concentration stack - raw data (uM)]");}
	else AppendToStack("Concentration stack - raw data (uM)", "FLIM Calcium concentration - raw data (in uM)");

	// Add scale bar and merge the calcium concentration
	selectWindow("FLIM Calcium concentration - raw data (in uM)");
	Border_size = floor(Border*height);
	//print("Border size: "+Border_size+" pixels");
	AddScaleBar("FLIM Calcium concentration - raw data (in uM)", Border_size);
	rename("FLIM Calcium concentration");
} //end of the if(Disp_CaImage == 1 || Disp_CaImage_merged == 1){


if(Disp_CaImage == 1){
	// Create and append the stack 
	if (j == 0) {
		selectWindow("FLIM Calcium concentration");
		run("Duplicate...", "title=[Concentration stack - RGB]");}
	else AppendToStack("Concentration stack - RGB", "FLIM Calcium concentration");
}


if(Disp_CaImage_merged == 1){
	selectWindow("FLIM total intensity");
	run("Duplicate...", "title=[FLIM total intensity - mask]");
	getMinAndMax(min, max);
	run("Divide...", "value=max");

	new_width = width + Border_size;
	run("Canvas Size...", "width=&new_width height=&height position=Center-Left zero");

	setColor(1);
	fillRect(width+1, 0, Border_size, height);
	MergeIntensityImage("FLIM Calcium concentration", "FLIM total intensity - mask");

	close("FLIM total intensity - mask");
	close("FLIM Up");
	close("FLIM Bottom");

	// Create and append the stack 
	if (j == 0) {
		selectWindow("FLIM Calcium concentration merged");
		run("Duplicate...", "title=[Concentration merged stack - RGB]");}
	else AppendToStack("Concentration merged stack - RGB", "FLIM Calcium concentration merged");	
}



// The raw data are now saved in the stacks so there is no need to keep it
close("FLIM lifetime image - raw data (in ns)")
;
close("FLIM Calcium concentration - raw data (in uM)");

// Close unnecessary images
if (Keep_nonResults == 0){
	close("FLIM stack");
	close("FLIM mask");}

// These files need to be renamed after the calculation in order to use token names for the calculation
if (Keep_nonResults == 1){
	selectWindow("FLIM stack");
	rename("FLIM stack - "+Short_filename);
	selectWindow("FLIM mask");
	rename("FLIM mask - "+Short_filename);}

if(Disp_IntImage == 0) close("FLIM total intensity");
if(Disp_IntImage == 1){
	selectWindow("FLIM total intensity");
	rename("CaliFLIM total intensity - "+Short_filename);
			// Create and append the lifetime stack 
	if (j == 0) {
		selectWindow("CaliFLIM total intensity - "+Short_filename);
		run("Duplicate...", "title=[Total intensity stack - raw data (ADC)]");}
	else AppendToStack("Total intensity stack - raw data (ADC)", "CaliFLIM total intensity - "+Short_filename);
}

if(Disp_TauImage == 1){
	selectWindow("FLIM lifetime image");
	rename("CaliFLIM lifetime image - "+Short_filename);}
	
if(Disp_TauImage_merged == 1){
	selectWindow("FLIM lifetime image merged");
	rename("CaliFLIM lifetime image merged - "+Short_filename);}

if(Disp_CaImage == 1){
	selectWindow("FLIM Calcium concentration");
	rename("CaliFLIM Calcium concentration - "+Short_filename);
}
	
if(Disp_CaImage_merged == 1){
	selectWindow("FLIM Calcium concentration merged");
	rename("CaliFLIM Calcium concentration merged - "+Short_filename);}


// Saving the images as selected by the user
if (Save_SingleImages_on == 1){
	//creates folder "output"

	if(Disp_TauImage == 1) 	SaveImage("CaliFLIM lifetime image - "+Short_filename, Output_dir_SingleImages);
	if(Disp_TauImage_merged == 1) SaveImage("CaliFLIM lifetime image merged - "+Short_filename, Output_dir_SingleImages);
	if(Disp_CaImage == 1) SaveImage("CaliFLIM Calcium concentration - "+Short_filename, Output_dir_SingleImages);

	if(Disp_CaImage_merged == 1) SaveImage("CaliFLIM Calcium concentration merged - "+Short_filename, Output_dir_SingleImages);
	if(Disp_IntImage == 1) SaveImage("CaliFLIM total intensity - "+Short_filename, Output_dir_SingleImages);
}


if (LargeBatchMode == false) setBatchMode("exit and display");

// Plot decays
if (Keep_nonResults == 1){
	Plot.create("FLIM decays - "+Short_filename, "Time (ns)", "Intensity (AU)");
	Plot.setLimits(0, t_window, -0.1, 1.1);
	if (IRFtype_choice != "None"){
		Plot.setColor("black");
		Plot.add("line", t_IRF, IRF_decay);
}
	Plot.setColor("red");
	Plot.add("line",t, Data_decay);
	Plot.show();}

if (LargeBatchMode == true) close("*FLIM*");
	
t_finish = getTime();
dt = (t_finish-t_start)/1000;
print("Time elapsed: "+dt+" s");

// End of for loop------------------------------------------------------------
}


print("-----------------------------------------------------------");

tf = getTime();
dt = (tf-t0)/1000; // in seconds
print("Total time elapsed: "+d2s(dt,1)+" seconds");

if (LargeBatchMode == true) setBatchMode("exit and display");
ArrangeWindows(Disp_IntImage, Disp_TauImage, Disp_TauImage_merged, Disp_CaImage, Disp_CaImage_merged);

// Save the log and results as txt files
if (Save_SingleImages_on == 1 || Save_Stacks_on == 1){
	selectWindow("Log");  //select Log-window
	saveAs("Text", Output_dir+"Log file");
	selectWindow("Results");  //select Log-window
	saveAs("Measurements", Output_dir+"Results.txt");}
	

// Display and save the stacks
// Lifetime data ----------------------------
if(Disp_TauImage == 1 || Disp_TauImage_merged == 1){
	selectWindow("Lifetime stack - raw data (ns)");
	setSlice(1);
	//setBatchMode("show");
	if(Save_Stacks_on == 1){
		selectWindow("Lifetime stack - raw data (ns)");
		saveAs("Tiff",Output_dir+"\\"+"Lifetime stack - raw data (ns)");}
}

if(Disp_TauImage == 1 && Save_Stacks_on == 1){
		selectWindow("Lifetime stack - RGB");
		saveAs("Tiff",Output_dir+"\\"+"Lifetime stack - RGB");}
	
if(Disp_TauImage_merged == 1 && Save_Stacks_on == 1){
	selectWindow("Lifetime merged stack - RGB");
	saveAs("Tiff",Output_dir+"\\"+"Lifetime merged stack - RGB");}

// Concentration data ----------------------------
if(Disp_CaImage == 1 || Disp_CaImage_merged == 1){
	selectWindow("Concentration stack - raw data (uM)");
	setSlice(1);
	//setBatchMode("show");
	if(Save_Stacks_on == 1){
		selectWindow("Concentration stack - raw data (uM)");
		saveAs("Tiff",Output_dir+"\\"+"Concentration stack - raw data (uM)");}
}

if(Disp_CaImage == 1 && Save_Stacks_on == 1){
		selectWindow("Concentration stack - RGB");
		saveAs("Tiff",Output_dir+"\\"+"Concentration stack - RGB");}
		
	
if(Disp_CaImage_merged == 1 && Save_Stacks_on == 1){
	selectWindow("Concentration merged stack - RGB");
	saveAs("Tiff",Output_dir+"\\"+"Concentration merged stack - RGB");}

// Totral intensity data --------------------------
if(Disp_IntImage == 1){
	selectWindow("Total intensity stack - raw data (ADC)");
	setSlice(1);
	//setBatchMode("show");
	if (Save_Stacks_on == 1){
		selectWindow("Total intensity stack - raw data (ADC)");
		saveAs("Tiff",Output_dir+"\\"+"Total intensity stack - raw data (ADC)");}
}



// -----------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------

// Functions


// ------------------------------------------------------------
// Function that opens the SDT file and read the appropriate Metadata
function OpenFLIMdata(filename){
	print(filename);

	if (endsWith(filename, ".sdt")){
		Ext.setSeries(0);
		Ext.setId(filename);
		Ext.getMetadataValue("bhfileHeader.noOfDataBlocks", n_DataBlock);
		n_DataBlock = parseInt(n_DataBlock);
		// Case where the IRF is saved as a multiblock file
		if (n_DataBlock > 1) {
			Block_number = newArray(n_DataBlock);
			for (s=0; s<n_DataBlock; s++) Block_number[s] = d2s(s, 0);

			SplitPath = split(filename, "\\");
			Dialog.create(GUI_name+" - Multi-block dataset");
			Dialog.addMessage("----- Choose data block to use -----");
			Dialog.addMessage(SplitPath[SplitPath.length-1]);
			Dialog.addChoice("Data block: ", Block_number, Block_number[n_DataBlock-1]);
			Dialog.show();
			Block_Choice = Dialog.getChoice();
			Ext.setSeries(parseInt(Block_Choice));}
		Ext.openImagePlus(filename);

		// Alternatively, exclude the possibility of opening multiple channel IRF
		//exit("Wrong format (please see manual to convert your IRF to a single channel IRF).");
		
	// Increment factor was not automatically corrected for in version of Bio-format before 5.0.1
	// So it previously needed to do it
	
	//Correct_IncrFactor = 0;
	//if (Correct_IncrFactor == 1) {
	//	Ext.getMetadataValue("MeasureInfo.incr", IncrFactor);
	//	IncrFactor = parseFloat(IncrFactor);
		//print("Increment factor: " +IncrFactor);

		// Correct for the increment factor
	//	run("Divide...", "value=&IncrFactor stack");}
		
		}
		
	//if (endsWith(filename, ".tif") || endsWith(filename, ".tiff")) {
	//	if (endsWith(filename, ".ome.tif")) {
	//		Ext.setSeries(0);
	//		Ext.setId(filename);
	//		Ext.openImagePlus(filename);}
	//	else open(filename);
	//	}

	if (endsWith(filename, ".tif") || endsWith(filename, ".tiff")) open(filename);

	run("32-bit");
	// This next lines are important when opening .tif files that have been saved in different ways
	getDimensions(width, height, channels, slices, frames);
	if (channels > frames) run("Re-order Hyperstack ...", "channels=[Frames (t)] slices=[Slices (z)] frames=[Channels (c)]");
	if (slices > frames) run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");

	// print("File opened successfully.");
}


// ------------------------------------------------------------
// GetTimeWindow gets the parameter t_window depending on the file type
function GetTimeWindow(filename){
		//print(filename);

	if (endsWith(filename, ".sdt")){
		Ext.setId(filename);		
		Ext.getMetadataValue("time base", t_window);
		t_window = parseFloat(t_window);}
		
	if (endsWith(filename, ".tif") || endsWith(filename, ".tiff")){
		open(filename);
		Filename_wo_path = getTitle;
		close();
		// Create the dialog window
		Dialog.create("Setting time window");
		Dialog.addMessage("----- Setting time window -----");
		Dialog.addMessage("File name: "+Filename_wo_path);
		Dialog.addNumber("Time window", 25, 2, 5, " ns");
		Dialog.show();
		t_window = Dialog.getNumber();}
	
return t_window;}



// ------------------------------------------------------------
// Calculate the total intensity and the mask
function CalculateMask(ImageName, Thresh_param){
	selectWindow(ImageName);
	//getDimensions(width, height, channels, slices, frames);
	//rename("Intensity image for masking");
	
	run("Duplicate...", "title=["+ImageName+" mask]");
	run("Gaussian Blur...", "sigma=1");
	run("Subtract Background...", "rolling=50");
	GetRidNegative(ImageName+" mask");
	run("Remove Outliers...", "radius=2 threshold=100 which=Bright");
	run("Remove Outliers...", "radius=2 threshold=100 which=Dark");

	setAutoThreshold("Li dark");
	setOption("BlackBackground", false);
	// Get the threshold defined by Li thresholding
	getThreshold(lower, upper);
	// Apply the thresholding parameter to the lower limit
	setThreshold(Thresh_param*lower, upper);

	// Convert to binary from the thresholding
	run("Make Binary");
	run("Erode");
	run("Dilate");
	run("Divide...", "value=255");
	run("Red");
	setMinAndMax(0, 0);
	
	//close("Intensity image for masking");
}

// ------------------------------------------------------------

// Function that may be applied if we decide to apply different masking for the dataset and the calibration files (Tau_min and Tau_max)
function CalculateMask2(ImageName){
	selectWindow(ImageName);
	//getDimensions(width, height, channels, slices, frames);
	//run("Z Project...", "start=1 stop="+frames+" projection=[Sum Slices]");
	//rename("Intensity image for masking");

	run("Duplicate...", "title=["+ImageName+" mask]");
	run("Gaussian Blur...", "sigma=2");
	setAutoThreshold("Default dark");
	//setAutoThreshold("Minimum dark");
	run("Convert to Mask");
	run("Divide...", "value=255");
	run("Red");
	setMinAndMax(0, 0);	

	//close("Intensity image for masking");
}

// ------------------------------------------------------------
// Calculate the average lifeitme from dataset
function CalculateAvLifetime(Image_stack, t_window, Remove_bg, BGbin_min, BGbin_max, Tau_irf, Apply_corr_factor, Apply_CMF3fusion) {

	selectWindow(Image_stack);
	getDimensions(width, height, channels, slices, frames);
	dt = t_window/frames; // dt is the same for all analysis windows T

	// Remove the background
	//if (Remove_background == 0) print("Background not removed.");
	if (Remove_bg == 1){
		RemoveBackgroundFromStack(Image_stack, BGbin_min, BGbin_max);
		close("Background_stack");
		close(Image_stack+" - Background_image");
		//print("Background removed (using bins between bin #"+BGbin_min+" and bin #"+BGbin_max+")");
	}

	// Calculate the TauCMM image for T=Tmax
	selectWindow(Image_stack);
	run("Duplicate...", "duplicate");
	rename("ImageStack T");
	CalculateTauCMM("ImageStack T", dt, Tau_irf);
	
	// Applying the correction factor
	if (Apply_corr_factor == true){
		T_corr = t_window - Tau_irf;
		n_iteration = 10;
		ApplyCMM_correction("ImageStack T - Lifetime image", T_corr, n_iteration);}
		
	GetRidNegative("ImageStack T - Lifetime image");
	close("ImageStack T");

	if (Apply_CMF3fusion == false){
		selectWindow("ImageStack T - Lifetime image");
		rename(Image_stack+" Lifetime image no mask");}

	if (Apply_CMF3fusion == true){
		// Calculate the TauCMM image for T=Tmax/2
		n_frames_T2 = frames/2;
		selectWindow(Image_stack);
		run("Duplicate...", "duplicate");
		rename("ImageStack T/2");
		run("Slice Remover", "first="+n_frames_T2+1+" last="+frames+" increment=1");
		CalculateTauCMM("ImageStack T/2", dt, Tau_irf);
	
		// Applying the correction factor
		if (Apply_corr_factor == true){
			T_corr = t_window/2 - Tau_irf;
			n_iteration = 10;
			ApplyCMM_correction("ImageStack T/2 - Lifetime image", T_corr, n_iteration);
		}
		
		GetRidNegative("ImageStack T/2 - Lifetime image");
		close("ImageStack T/2");
	
		// Calculate the TauCMM image for T=Tmax/4
		n_frames_T3 = frames/4;
		selectWindow(Image_stack);
		run("Duplicate...", "duplicate");
		rename("ImageStack T/4");
		run("Slice Remover", "first="+n_frames_T3+1+" last="+frames+" increment=1");
		CalculateTauCMM("ImageStack T/4", dt, Tau_irf);
	
		// Applying the correction factor
		if (Apply_corr_factor == true){
			T_corr = t_window/4 - Tau_irf;
			n_iteration = 10;
			ApplyCMM_correction("ImageStack T/4 - Lifetime image", T_corr, n_iteration);
		}
		
		GetRidNegative("ImageStack T/4 - Lifetime image");
		close("ImageStack T/4");
	
		// calculate the weighting factors and carry out fusion
		T1 = t_window - Tau_irf;
		T2 = t_window/2 - Tau_irf;
		T3 = t_window/4 - Tau_irf;
	
		// Using T/2 lifetime image as seed 
		CalculateWeightingImage("ImageStack T/2 - Lifetime image", T1, T2);
		rename("W12");
	
		CalculateWeightingImage("ImageStack T/2 - Lifetime image", T2, T3);
		rename("W23");
	
		imageCalculator("Multiply create 32-bit", "ImageStack T - Lifetime image","W12");
		rename("WTau1");
	
		imageCalculator("Subtract create 32-bit", "W23","W12");
		rename("WTau2");
		imageCalculator("Multiply 32-bit", "WTau2", "ImageStack T/2 - Lifetime image");
	
		selectWindow("W23");
		run("Duplicate...", "title=WTau3");
		run("Subtract...", "value=1");
		run("Multiply...", "value=-1");
		imageCalculator("Multiply 32-bit", "WTau3", "ImageStack T/4 - Lifetime image");
	
		// Add them all up
		imageCalculator("Add 32-bit", "WTau1", "WTau2");
		imageCalculator("Add 32-bit", "WTau1", "WTau3");
	
		selectWindow("WTau1");
		rename(Image_stack+" Lifetime image no mask");
	
		close("W12");
		close("W23");
		close("WTau2");
		close("WTau3");
	
		close("ImageStack T - Lifetime image");
		close("ImageStack T/2 - Lifetime image");
		close("ImageStack T/4 - Lifetime image");
	}

}

// ------------------------------------------------------------
// Calculate the weighting image from the TauImage seed using the Tau_cutoff formulae
function CalculateWeightingImage(LifetimeImage, T1, T2){

	//A = 2.832;
	//C = 0.02309;
	
	A = 3.218;
	C = 0.07339;
	Tau_cutoff = sqrt(T1*T2*C/A);
	//print("Calculating the Weightng image for CMF3 using Tau_cutoff= "+Tau_cutoff+ " ns");
	
	selectWindow(LifetimeImage);
	run("Duplicate...", "title=WeightingImage");
	
	run("Subtract...", "value="+Tau_cutoff);
	run("Divide...", "value="+Tau_cutoff);
	run("Multiply...", "value=-20");
	run("Exp");
	run("Add...", "value=1");
	run("Reciprocal");

}


// ------------------------------------------------------------
// TauCMM lifetime image
function CalculateTauCMM(Image_stack_name, dt, Tau_irf){
	
	// Calculate lifetime image
	selectWindow(Image_stack_name);
	
	getDimensions(width, height, channels, slices, frames);
	//print("Calculating Tau_CMM on "+frames+" frames.");
	
	newImage("FLIM t-stack", "32-bit black", width, height, frames);
	
	// Creates the t image stack varying from 0 to 25 ns - dt (256 bins)
	for (i=0;i<frames;i++){
		setSlice(i+1); 
		t = dt*i;
		run("Add...", "value=t slice");
	}
	
	// Calculate the STd image
	imageCalculator("Multiply create stack", Image_stack_name,"FLIM t-stack");
	rename("FLIM STd image stack");
	run("Z Project...", "start=1 stop="+frames+" projection=[Sum Slices]");
	rename("FLIM STd image");
	
	// Calculate the Sd image
	selectWindow(Image_stack_name);
	run("Z Project...", "start=1 stop="+frames+" projection=[Sum Slices]");
	rename("FLIM Sd image");
	
	// Divide the 2 to obtain the lifetime image
	imageCalculator("Divide create 32-bit", "FLIM STd image","FLIM Sd image");
	rename(Image_stack_name+" - Lifetime image");
	
	// Half-bin correction
	selectWindow(Image_stack_name+" - Lifetime image");
	Half_bin = dt/2;
	run("Add...", "value="+Half_bin);
	
	// Remove IRF lifetime
	selectWindow(Image_stack_name+" - Lifetime image");
	run("Subtract...", "value="+Tau_irf);
	
	close("FLIM STd image stack");
	close("FLIM STd image");
	close("FLIM t-stack");
	close("FLIM Sd image");

}

// ------------------------------------------------------------
// CMM correction based on iterative calculation
function ApplyCMM_correction(Filename, T_corr, n_iteration){

	//t1 = getTime();
	selectWindow(Filename);
	run("Duplicate...", " ");
	rename("Pre-correction lifetime image");
	
	for (i=0;i<n_iteration;i++){
		selectWindow(Filename);
		run("Duplicate...", "title=ExpImage");
		close(Filename);
		
		run("Reciprocal");
		run("Multiply...", "value="+T_corr);
		run("Multiply...", "value=-1");
		run("Exp");
	
		run("Duplicate...", "title=NumeratorImage");
		run("Multiply...", "value="+T_corr);
	
		selectWindow("ExpImage");
		rename("DenominatorImage");
		run("Subtract...", "value=1");
		run("Multiply...", "value=-1");
		imageCalculator("Divide create 32-bit", "NumeratorImage","DenominatorImage");
		rename("Lifetime correction image");
	
		imageCalculator("Add create 32-bit", "Pre-correction lifetime image", "Lifetime correction image");
		rename(Filename);
	
		close("Lifetime correction image");
		close("DenominatorImage");
		close("NumeratorImage");}
	
	close("Pre-correction lifetime image");
	//t2 = getTime();
	//Duration = (t2-t1)/1000;
	//print("Time elapsed applying correction: "+Duration+" s");

}




// ------------------------------------------------------------
// Function to add the scale bar on the right side of the image
function AddScaleBar(Image_name, Border_size) {

	selectWindow(Image_name);
	getDimensions(width, height, channels, slices, frames);
	new_width = width + Border_size;

	// Empirically tried that and it worked
	zoom = 1.6*height/256;

	run("Duplicate...", "title=[Bigger image]");
	run("Canvas Size...", "width=&new_width height=&height position=Center-Left zero");
	run("Calibration Bar...", "location=[Upper Right] fill=Black label=White number=5 decimal=1 font=10 zoom=&zoom");
	rename(Image_name+" with scale bar");
	close("Bigger image");
}



// ------------------------------------------------------------
// Function that merges the colored image with the intensity image
function MergeIntensityImage(Image_name, Intensity_Image) {
	selectWindow(Image_name);
	run("Duplicate...", "title=["+Image_name+" merged]");
	run("RGB Color");
	run("Split Channels");

	imageCalculator("Multiply create 32-bit", Image_name+" merged (red)", Intensity_Image);
	rename("Merged red");
	run("8-bit");

	imageCalculator("Multiply create 32-bit", Image_name+" merged (green)", Intensity_Image);
	rename("Merged green");
	run("8-bit");

	imageCalculator("Multiply create 32-bit", Image_name+" merged (blue)", Intensity_Image);
	rename("Merged blue");
	run("8-bit");

	run("Merge Channels...", "c1=[Merged red] c2=[Merged green] c3=[Merged blue] create");
	rename("Merged 8-bits");

	run("RGB Color");
	rename(Image_name+" merged");

	close("Merged red");
	close("Merged green");
	close("Merged blue");
	close(Image_name+" merged (red)");
	close(Image_name+" merged (green)");
	close(Image_name+" merged (blue)");
	close("Merged 8-bits");
}



// ------------------------------------------------------------
// Function that gets rid of the NaN values in the image by 0
function GetRidNegative(ImageName) {
	newValue = 0;
	counter = 0;

	selectWindow(ImageName);
	for (y = 0; y < getHeight(); y++){
        	for (x = 0; x < getWidth(); x++){
                	p = getPixel(x,y);
                	if (p<0) {
                        	setPixel(x, y, newValue);
                        	counter++;
                	}
        	}
	}

//print(Image+" - Get rid of negative values");
//print("" + counter + " pixels replaced"); 
}


// ------------------------------------------------------------
// Save the image in the right folder and name and closes it
function SaveImage(Image_name, path){
	selectWindow(Image_name);
	run("Duplicate...", "title=[Save image]");
	SaveFilename = path+"/"+Image_name+".png";
	saveAs("png", SaveFilename);
	close();
	print("Saved:");
	print(SaveFilename);
}	


// ------------------------------------------------------------
// Function that makes the list of files containing the token name Data_name
function Make_data_file_list(dir,Data_name, Extension){
	//print("Making list of datasets");
	AllFile_list = getFileList(dir);
	//Array.show(AllFile_list);
	DataFile_list = newArray(0);

	for (i=0;i<AllFile_list.length;i++) {
		sdt_file = endsWith(AllFile_list[i], Extension);
		contains_name = indexOf(AllFile_list[i], Data_name)>0 || startsWith(AllFile_list[i], Data_name);
		if (sdt_file && contains_name == true) DataFile_list = Array.concat(DataFile_list, dir + AllFile_list[i]);
	}
	
	// Exit if the file is not found

	if (lengthOf(DataFile_list) == 0) exit("No files containing *"+Data_name+"*"+Extension+" were found.");
	
	return DataFile_list;
}


// ------------------------------------------------------------
// Remove the background from a stack
function RemoveBackgroundFromStack(Image_stack, BGbin_min, BGbin_max){
	selectWindow(Image_stack);
	run("32-bit");
	rename(Image_stack+"_original");
	run("Make Substack...", "slices="+BGbin_min+"-"+BGbin_max);
	rename("Background_stack");
	run("32-bit");
	
	run("Z Project...", "start="+BGbin_min+" stop="+BGbin_max+" projection=[Average Intensity]");
	rename(Image_stack+" - Background_image");
	
	run("Gaussian Blur...", "sigma=1");
	run("Median...", "radius=5");
	imageCalculator("Subtract create 32-bit stack", Image_stack+"_original",Image_stack+" - Background_image");
	
	rename(Image_stack);
	
	close(Image_stack+"_original");	
	close("Background_stack");
	close(Image_stack+" - Background_image");
	
	//print("Background removed (using bins between bin #"+BGbin_min+" and bin #"+BGbin_max+")");
}



// ------------------------------------------------------------
// Get the decay data from the image stack
function GetDecayData(Stack_name, Remove_bg, BGbin_min, BGbin_max){
	selectWindow(Stack_name);
	getDimensions(width, height, channels, slices, frames);
	run("Duplicate...", "title=[Sum_stack] duplicate frames=1-frames");
	run("32-bit");

	if (Remove_bg == 1) {
		//print("GetDecayData: Removing background");
		RemoveBackgroundFromStack("Sum_stack", BGbin_min, BGbin_max);
		//print("Done");
		}

	
	run("Bin...", "x=&width y=&height z=1 bin=Sum");

	// Get the profile in an array
	Z_profile = newArray(frames);
	for (i = 1; i<frames+1; i++){
		setSlice(i);
		Z_profile[i-1] = getPixel(0,0);}

	// Normalise the data
	Array.getStatistics(Z_profile, min, max, mean, stdDev);
	for (i=0;i<Z_profile.length;i++) Z_profile[i] = Z_profile[i]/max;
	close("Sum_stack");

	return Z_profile;
}


// ------------------------------------------------------------
// Get the temporal position @ a certain threshold for shift estimation
function GetTimeThresh(Time,Decay,Thresh){

	i = 1;
	while (Decay[i] < Thresh) i++;
	
	// print(Decay[i]);
	// print(i);
	// print(Time[i]);
	// print(Decay[i-1]);
	
	// Calculate T_thresh
	if (Decay[i] == Thresh) {
		T_thresh = Time[i];}
	else {
	// Linear interpolation
	a = (Decay[i] - Decay[i-1])/(Time[i] - Time[i-1]);
    b = Decay[i] - a*Time[i];
    T_thresh = (Thresh - b)/a;}

	return T_thresh;
}


// ------------------------------------------------------------
// Calculate the SHift correction based on the IRF standard deviation (empirical...)
function Shift_correction_factor(t_IRF, IRF_decay){
	print("Gaussian fitting (no offset)...");
	Fit.doFit("Gaussian (no offset)", t_IRF, IRF_decay);
	print("A = "+d2s(Fit.p(0),2)+" (AU)");
	print("t0 = "+d2s(Fit.p(1),3)+" ns");
	print("Sigma = "+d2s(Fit.p(2),3)+" ns");

	// Empirical equation
	Shift_correction = 0.4*Fit.p(2)+0.02;
	print("Shift correction: "+d2s(Shift_correction,3)+" ns");
	  
	return Shift_correction;
}


// ------------------------------------------------------------
// Calculate the bin position of the bottom of the rising edge from T_thresh, this can be used to determine the bins to use for background estimation
function TimeThresh2Bin_position(T_thresh, t_window, frames){
	dt = t_window/frames;
	Bin_position = floor(T_thresh/dt);
	if (Bin_position > 2) Bin_position = Bin_position - 2;
	else Bin_position = 0;
	print("Bin position of Edge threshold: "+Bin_position);

	return Bin_position;
}

// ------------------------------------------------------------
// This function gets the position of the background estimation by calculating the derivative of the curve and using a threshold
function GetBackgroundPosition(Decay){
	DecayLength = lengthOf(Decay);
	Y_1 = Array.slice(Decay,0,DecayLength - 2);
	Y_2 = Array.slice(Decay,2,DecayLength);

	dY = newArray(DecayLength-2);
	for (i=0;i<DecayLength-2;i++) dY[i] = Y_2[i] - Y_1[i]; // derivative so that it is around zero

	X = Array.getSequence(DecayLength - 2); 
	Position = GetTimeThresh(X,dY,0.02);
	Position = floor(0.9*Position) - 2;
	print("Background bins from 1 to "+Position);
	
	return Position;
}

// ------------------------------------------------------------
// This function appends an image at the end of a stack
function AppendToStack(StackName, ImageName) {
	selectWindow(ImageName);
	run("Select All");
	run("Copy");
	run("Select None");
	selectWindow(StackName);
	getDimensions(width, height, channels, slices, frames);
	setSlice(slices);
	run("Add Slice");
	run("Paste");
	run("Select None");
	setBatchMode("hide");
}


// ------------------------------------------------------------
// This function arranges the windows in a sensible way
function ArrangeWindows(Disp_IntImage, Disp_TauImage, Disp_TauImage_merged, Disp_CaImage, Disp_CaImage_merged){
	SC_width = screenWidth;
	SC_height = screenHeight;

	//print(SC_width);
	//print(SC_height);

	x0 = floor(SC_width/10);
	y0 = floor(SC_height/8);
	width0 = 0;
	height0 = 0;

	x1 = x0;
	y1 = y0;
	width1 = 0;
	height1 = 0;

	x2 = x0;
	y2 = y0;
	width2 = 0;
	height2 = 0;
	

	if(Disp_IntImage == 1){
		selectWindow("Total intensity stack - raw data (ADC)");
		setLocation(x0, y0);
		getLocationAndSize(x0, y0, width0, height0);
		x1 = x0 + width0;
		y1 = y0;
		x2 = x1;
		y2 = y1;}

	if(Disp_TauImage == 1 || Disp_TauImage_merged == 1){
		selectWindow("Lifetime stack - raw data (ns)");
		setLocation(x0 + width0, y0);
		getLocationAndSize(x1, y1, width1, height1);
		x2 = x1 + width1;
		y2 = y1 + height1;}

	if(Disp_TauImage == 1){
		selectWindow("Lifetime stack - RGB");
		setLocation(x1 + width1, y0);
		getLocationAndSize(x2, y2, width2, height2);}

	if(Disp_TauImage_merged == 1){
		selectWindow("Lifetime merged stack - RGB");
		setLocation(x2 + width2, y0);}
	

	if(Disp_CaImage == 1 || Disp_CaImage_merged == 1){
		selectWindow("Concentration stack - raw data (uM)");
		setLocation(x1, y1 + height1);
		getLocationAndSize(x1, y1, width1, height1);
		x2 = x1 + width1;
		y2 = y1 + height1;}

	if(Disp_CaImage == 1){
		selectWindow("Concentration stack - RGB");
		setLocation(x1 + width1, y1);
		getLocationAndSize(x2, y2, width2, height2);
		x2 = x1 + width1;
		y2 = y1 + height1;}

	if(Disp_CaImage_merged == 1){
		selectWindow("Concentration merged stack - RGB");
		setLocation(x2 + width2, y1);}
}


// ------------------------------------------------------------
// This function finds the position of the maximum in the decay

function FindPositionOfMaximum(t, Decay) {

	RankedPositions = Array.rankPositions(Decay);
	RankedPositions = Array.reverse(RankedPositions);
	t_max = t[RankedPositions[0]];
	
	return t_max	
}


// Function that excludes the IRF datasets from the list of data to analyse
//function ExcludeIRF_dataset(AllFile_list, IRF_filename, IRF_Calib_filename){
//
//	DataFile_list = newArray(0);
//	for (i=0;i<AllFile_list.length;i++){
//		if (AllFile_list[i] != IRF_filename && AllFile_list[i] != IRF_Calib_filename){
//			DataFile_list = Array.concat(DataFile_list, AllFile_list[i]);}
//	}
	
//return DataFile_list;}


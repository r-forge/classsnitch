
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="http://<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="http://r-forge.r-project.org/"><img src="http://<?php echo $themeroot; ?>/imagesrf/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->
<!— Header —>
<div class="header">
<h1 align="center"><font face="verdana">classSNitch</font></h1>
<h2 align="center"><font face="Geneva">detection and classification of RNA structure change in chemical mapping data</font></h2>
</div>
<hr>

<!—  RiboSNitch Image —>
<table><tr>
<td><figure>
<img id="im1" src="images/riboSNitches.png" width="390" height="417" alt="RiboSNitches"/>
</figure></td>
<td>
<h3 align="center">About</h3>
<p  align="justify">Mutations in RNA will create a riboSNitch, if important structural elements are disrupted. Recent ultra-high throughput techniques, such as SHAPE-MaP and PARS, enable the collection of structural RNA information on a genome-wide scale. With the ability to gather genome-wide structural information on RNA, it is important to accurately classify these structural data in order to identify those structural changes that result in a phenotypic outcome. We have developed an automated approach to classify structure change in SHAPE data. This method utilizes random forest classification on a set of pattern and magnitude parameters from the mutate and map SHAPE data set (or another user specified data set) to build a classifier. The classifier is then used to identify structure change in other SHAPE traces. Enabling scientists to identify structure change may help guide experiments that examine RNA structure and its role in biological processes.</p>
</td>
</tr></table>
<hr>

<!—  Package Info —>
<table>
<tr valign="top">
<td width="25%"><figure>
<h3>Authors</h3>
</figure></td>
<td  width="75%">
<br>
Chanin Tolson and Alain Laederach
</td>
</tr >
<tr valign="top">
<td width="25%"><figure>
<h3>Download</h3>
</figure></td>
<td  width="75%">
<br>
Package source (use this for Linux / Solaris): <a href="docs/classSNitch_0.13.0.tar.gz" download>classSNitch_0.13.0.tar.gz</a><br>
Windows binary: <a href="docs/classSNitch_0.13.0.zip" download>classSNitch_0.13.0.zip</a><br>
Repository: classSNitch is also available on <a href="http://r-forge.r-project.org/projects/classsnitch/">R-Forge</a>
</td>
</tr>
<tr valign="top">
<td width="25%"><figure>
<h3>Installation</h3>
</figure></td>
<td  width="75%">
<br>
Linux: At the command line type "R CMD INSTALL classSNitch_0.13.0.tar.gz".<br><br>
Windows: From the pull-down menu, click on "Packages->Install Packages from local zip file", and then select the downloaded file docs/classSNitch_0.13.0.zip. <br><br>
R install: install.packages("classSNitch", repos="http://R-Forge.R-project.org", type="source")<br><br>
Make sure R is installed on your computer. It is freely available on <a href="http://www.r-project.org."></a>R-Forge<br>
You will also need randomForest, dtw, ROCR, and gplots packages (available from the <a href="https://cran.r-project.org/">Comprehensive R Archive Network)</a>.

</td>
</tr>
<tr valign="top">
<td width="25%"><figure>
<h3>Getting Started</h3>
</figure></td>
<td  width="75%">
<br>
Load classSNitch package: library("classSNitch")<br>
Get list of help pages: help(package="classSNitch")<br><br>
Documentation: <a href="docs/classSNitch.pdf">Reference Manual</a>
</td>
</tr>
<tr valign="top">
<td width="25%"><figure>
<h3>Examples</h3>
</figure>
</td>
<td  width="75%">
<br>
#get change parameters<br>
data("shape_ex")<br>
sample_shape = shape_ex<br>
sample = getChangeParams(sample_shape[2:nrow(sample_shape),], base=sample_shape[1,], trim=5)<br>
<br>
#predict change<br>
cr = classifyRNA(classes=2)<br>
cr_pred = predict(cr, sample, type="response")
</td>
</tr>
<tr valign="top">
<td width="25%"><figure>
<h3>Contact</h3>
</figure></td>
<td  width="75%">
<br>
Send an email to <a href="mailto:ctolson@email.unc.edu">Chanin Tolson</a>. We welcome any questions, comments or suggestions. If you have any publications using classSNitch, we would be interested in seeing how the package is being used.
</td>
</tr>
</table>

</body>
</html>

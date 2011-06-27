
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='http://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="/"><img src="<?php echo $themeroot; ?>/images/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>

<!-- end of project description -->

<p> [27 June 2011] Note that the Repitools projects is now being hosted by Bioconductor, from version <a href="http://www.bioconductor.org/packages/2.9/bioc/html/Repitools.html">2.9 onwards</a>.</p><br>

<FONT COLOR="gray">

<p> User's guide in PDF format can be downloaded <a href="RepitoolsManual.pdf">here</a>. </p>

<p> 'RepitoolsExamples' package containing example data can be downloaded from <a href="http://129.94.136.7/file_dump/mark/RepitoolsExamples_1.0.10.tar.gz">here</a>. </p>

<p> Users can join the mailing list <a href="https://lists.r-forge.r-project.org/cgi-bin/mailman/listinfo/repitools-help">here</a>.

<p> The R package can be downloaded from <a href="http://r-forge.r-project.org/R/?group_id=716">here</a> or installed directly using: <br> <br>
<code>install.packages("Repitools",repos="http://r-forge.r-project.org", type="source")</code></p>

<p> The <strong>project summary page</strong> you can find <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. </p>

</body>
</html>

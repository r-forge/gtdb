
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

<h2>Welcome to the Genotype Database Framework, GT.DB</h2>

<style type="text/css">p.big {line-height: 1.5}</style>
<p class="big">
GT.DB is a database framework for efficiently managing genome-wide
genotype datasets (i.e. 10<sup>6</sup> or more biallelic markers, 1000s of
samples), as well as underlying raw data and complex phenotype data,
tightly integrated with R for analysis.  It is based on components
developed at Perlegen Sciences, where a similar framework was used to
manage nearly 10<sup>11</sup> individual genotypes with underlying data,
representing more than 100,000 samples and more than 100 experimental
datasets.
</p>

<!-- end of project description -->

<p> The <strong>project summary page</strong> can be found <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. </p>

<h3>Getting Started</h3>

<ul>
  <li>Install the gt.db package
<pre>
install.packages("gt.db", repos="http://R-Forge.R-project.org")
library(gt.db)
</pre>
  </li>
  <li>Create a temporary SQLite database
<pre>
library(SQLite)
fname <- tempfile()
dbh <- dbConnect(dbDriver('SQLite'), fname)
</pre>
  </li>
  <li>Create the gt.db tables, and load up the demo datasets
<pre>
use.gt.db(dbh)
init.gt.db()
.gt.db.options(db.mode='hex')
demo('setup.gt.demo')
</pre>
  </li>
  <li>Explore the code examples!
<pre>
example(prcomp.gt.data)
example(xyplot.gt.data)
</pre>
  </li>
</ul>

</body>
</html>

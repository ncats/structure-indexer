Structure Indexer
=================

This is a self-contained structure indexer that uses Lucene as the
underlying storage and indexing engine. The indexing architecture
is based on an inverted indexing technique developed within NCATS.
Please consult the wiki (coming soon) for additional details.

Building
========

If you're building for the first time, you'll need to add the ```jchem.jar```
to your local maven repository. To do this, execute the following
command:

```
mvn install:install-file \
  -Dfile=lib/jchem.jar \
  -DgroupId=chemaxon \
  -DartifactId=jchem \
  -Dversion=3.2.12 \
  -Dpackaging=jar \
  -DgeneratePom=true
```

Then simply do ```mvn package``` to build the ```structure-indexer```
jar file. The ```bin``` directory contains the following wrapper scripts:

```indexer``` is the main driver that is used to build the index.
See ```indexer -h``` for complete usage. Here is an running example:

```
indexer index_dir BindingDB2D.sdf
```

```searcher``` is the client driver that provides a command-line interface
for searching and filtering. See ```searcher -h``` for complete usage. For
example, consider the following command:

```
searcher  -fsmiles -F_natoms=20:22 -F_molwt=280.:300. -F_source=BindingDB2D -s sim -t.9 idx "N1c2ccccc2NC(=O)c2cccnc12"
```

This example performs similarity searching against the index ```idx```
for the given structure with a Tanimoto cutoff of 0.9, number of atoms
in the range [20,22], molecular weight in the range [280, 300], 
only from the source ```BindingDB2D```, and outputs the matches as SMILES
format.

Note that there is currently a security alert https://github.com/ncats/structure-indexer/network/alert/pom.xml/org.apache.lucene:lucene-core/open but we do not expect a problem as we’re not using the feature that has the vulnerability

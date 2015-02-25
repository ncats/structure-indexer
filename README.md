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
jar file. The script ```indexer``` provides a simple wrapper to generate
a simple test search based on given molecule files; e.g., 

```
./indexer index_dir BindingDB2D.sdf
```

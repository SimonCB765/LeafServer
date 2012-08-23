from django.db.models.signals import post_syncdb
from django.db import connection, transaction

import Leaf.models

# Set name field to be BINARY, to force case-sensitive comparisons
def set_name_to_binary(sender, **kwargs):
    cursor = connection.cursor()
    cursor.execute('ALTER TABLE Leaf_proteininformation MODIFY chain varchar(5) BINARY NOT NULL')
    cursor.execute('ALTER TABLE Leaf_chaintype MODIFY chain varchar(5) BINARY NOT NULL')
    cursor.execute('ALTER TABLE Leaf_similarity MODIFY chainA varchar(5) BINARY NOT NULL')
    cursor.execute('ALTER TABLE Leaf_similarity MODIFY chainB varchar(5) BINARY NOT NULL')
    cursor.execute('ALTER TABLE Leaf_representative MODIFY nonreprChain varchar(5) BINARY NOT NULL')
    cursor.execute('ALTER TABLE Leaf_representative MODIFY reprChain varchar(5) BINARY NOT NULL')
    cursor.execute('ALTER TABLE Leaf_entryrepresentative MODIFY nonreprChain varchar(5) BINARY NOT NULL')
    cursor.execute('ALTER TABLE Leaf_entryrepresentative MODIFY reprChain varchar(5) BINARY NOT NULL')
    transaction.commit_unless_managed()

post_syncdb.connect(set_name_to_binary, sender=Leaf.models)
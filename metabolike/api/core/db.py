from .config import db_conf
from metabolike.db.metacyc import MetaDB

db = MetaDB(
    db_conf.metabolike_db_uri,
    db_conf.metabolike_db_user,
    db_conf.metabolike_db_password.get_secret_value(),
)
db.use_database(db_conf.metabolike_db_name)

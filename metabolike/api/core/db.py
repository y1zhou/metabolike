from metabolike.db.metacyc import MetaDB

from .config import db_conf

db = MetaDB(
    db_conf.metabolike_db_uri,
    db_conf.metabolike_db_user,
    db_conf.metabolike_db_password.get_secret_value(),
    db_conf.metabolike_db_name,
)

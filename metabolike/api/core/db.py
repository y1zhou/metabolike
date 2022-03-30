from metabolike.db.metacyc import MetacycClient

from .config import db_conf

db = MetacycClient(
    db_conf.metabolike_db_uri,
    db_conf.metabolike_db_user,
    db_conf.metabolike_db_password.get_secret_value(),
    db_conf.metabolike_db_name,
)
